"""
This module contain the alignment functions for reads to haplotype,
genotyping or realignment process
"""
import os
import ctypes

import common as com # Same common functions
import datum  as DM   # The global common datum
from read import Read

# DM.CommonDatum().hashmer is defualt to be 7
COMDM  = DM.CommonDatum()
ALIMER = 15 # Fix by the alignment algorithm in 'align.so'
# The diretory of this program
dir    = os.path.dirname(os.path.abspath(__file__))
align  = ctypes.CDLL(dir + '/align.so')  # The alignment module written by C

# Create a class to get the alignment information which record in a 
# c structure data type in align.fastAlignmentRoutine
# Thinks for the help of ``http://blog.csdn.net/joeblackzqq/article/details/10441017``
class AlignTagPointer(ctypes.Structure):
	_fields_ = [("score", ctypes.c_int), ("pos", ctypes.c_int)]
align.fastAlignmentRoutine.restype = ctypes.POINTER(AlignTagPointer)


def alignReadToHaplotype(haplotype, reads_collection, bam_stream):
    """ Learn from Platypus.

    This is the most basic and important part of AsmVar2. This function decides
    where to anchor the read to the specified haplotype, and calls the banded-
    alignment function -- fastAlignmentRoutine function in 'align.c'. 

    Args:
        `haplotype`:
        `reads_collection`: A hash to reads. It's a contianer of reads, use it 
                            to prevent recalculting the hash-sequence for the 
                            same reads. This could save the running time
    Return:

    Requires pysam
    """

    if haplotype.seq_hash is None:
        haplotype.seq_hash = com.SeqHashTable(haplotype.sequence)

    # 0-index system
    start, end = haplotype.hap_start - 1, haplotype.hap_end
    for r in bam_stream.fetch(haplotype.chrom, start, end):

        # The uniq-hash Id for identifying each read by the name and seq
        r_identify = hash((r.qname, r.query))
        if r_identify not in reads_collection:
            reads_collection[r_identify] = Read(r)

        # The alignemnt position in pysam is 0-base, shift r.pos to be 1-base
        score = singleRead2Haplotype(haplotype, 
                                     reads_collection[r_identify],
                                     r.pos + 1) # read's alignment position

    
def singleRead2Haplotype(haplotype, read, read_align_pos, do_calcu_flank_score = None):
    """
    Mapping a single read to a specified haplotype sequence. 

    This is used to find the best anchor position in the haplotype sequence
    for the specified read.

    Args:
        'haplotype': value of ``Haplotype`` class
        'read'     : value of ``Read`` class, which contain read's information 
                     collecting from ``AlignedRead`` of pysam by reading bam
                     files and a hash table for the read sequence

    Return: The best alignement score
    """
    # hash the halotype.sequence and read.seqs. And just do it here, and 
    # just do it once!!!
    if read.hash is None:
         read.seq_hash = com.SeqHashTable(read.seqs, COMDM.hashmer)
    if haplotype.seq_hash is None:
        haplotype.seq_hash = com.SeqHashTable(haplotype.sequence, COMDM.hashmer)
        haplotype.gap_open = com.set_gap_open_penalty(haplotype.sequence, 
                                                      COMDM.homopol_penalty)

    # Firstly, we use seq_hash to find the most possible anchor position
    max_hit = 0
    idx     = {}
    # Scan read sequence hash table
    for i, id in enumerate(read.seq_hash.hash_pointer): 

        if id in haplotype.seq_hash.hash_table:
            # Index in list of haplotype.seq_hash[id]
            # [BUG] This could be a bug if we hit tandem repeat regions, we may
            # always get the first position at the begin, because we don't use 
            # any strategies to select the most likely start from all the index 
            # in this list(haplotype.seq_hash.hash_table[id])!
            idx[id] = idx.get(id, -1) + 1     # Start from index 0
            pos_idx = haplotype.seq_hash.hash_table[id][idx[id]]
            haplotype.map_depth[pos_idx] += 1 # Record the mapping depth
            if haplotype.map_depth[pos_idx] > max_hit:
                max_hit = haplotype.map_depth[pos_idx]

    # max_hit > 0 means we can find some positions that read could anchor in 
    # haplotype by hash searching.
    
    hap_len_for_align = len(read) + ALIMER
    aln1 = ['\0' for i in range(2 * len(read) + ALIMER)]
    aln2 = ['\0' for i in range(2 * len(read) + ALIMER)]
    best_ali_score = 1000000 # A big enough bad score
    best_ali_pos   = -1
    firstpos       = 0
    if max_hit > 0:
        # Go through all the positions of haplotype.sequence 
        for i, d in enumerate(haplotype.map_depth):
            # Now let's find the best anchor position!
            # 'i' could represent the index of read in haplotype
            is_inside = i + hap_len_for_align <= len(haplotype)
            if d == max_hit and is_inside: 
                # Just start from the most possible position
                read_start_in_hap = max(0, i - 8) # 0-base system
                s = read_start_in_hap
                e = s + hap_len_for_align
                ali = align.fastAlignmentRoutine(haplotype.sequence[s:e], 
                                                 read.seqs, 
                                                 read.qual, 
                                                 hap_len_for_align, 
                                                 len(read), 
                                                 COMDM.gap_extend,
                                                 COMDM.nucprior,
                                                 haplotype.gap_open[s:e],
                                                 aln1, aln2)
                # ali.contents.pos == -1 means None in ``fastAlignmentRoutine``!
                # so that we should not refresh ``firstpos``
                if ali.contents.pos != -1: 
                    firstpos = ali.contents.pos

                # calculate contribution to alignment score of mismatches
                # and indels in flank, and adjust score. short circuit if
                # calculation is unnecessary
                if (do_calcu_flank_score and
                   haplotype.buffer_size and
                   ali.contents.score > 0):
                    ali.contents.score -= align.calculateFlankScore(
                        len(haplotype), 
                        haplotype.buffer_size,
                        read.qual,
                        haplotype.gap_open,
                        COMDM.gap_extend,
                        COMDM.nucprior,
                        read_start_in_hap + firstpos,
                        aln1, aln2)

                if ali.contents.score < best_ali_score:
                    best_ali_score = ali.contents.score
                    best_ali_pos   = i 

                    # Short-circuit this loop if we find an exact match
                    # (0 is the best score, means exact match)
                    if best_ali_score == 0: 
                        return best_ali_score

    # Try original mapping position. If the read is past the end of the 
    # haplotype then don't allow the algorithm to align off the end of 
    # the haplotype. This will only happen in the case of very large deletions.
    read_start_in_hap = min(len(haplotype) - hap_len_for_align, 
                            read_align_pos - haplotype.hap_start)

    if read_start_in_hap != best_ali_pos:
        read_start_in_hap = max(0, i - 8)
        s = read_start_in_hap 
        e = s + hap_len_for_align 
        ali = align.fastAlignmentRoutine(haplotype.sequence[s:e], 
                                         read.seqs, read.qual, 
                                         hap_len_for_align, len(read),
                                         COMDM.gap_extend, 
                                         COMDM.nucprior, 
                                         haplotype.gap_open[s:e],
                                         aln1, aln2)

        if ali.contents.pos != -1:
            firstpos = ali.contents.pos

        if (do_calcu_flank_score and haplotype.buffer_size and 
           ali.contents.score > 0): 
            ali.contents.score -= align.calculateFlankScore(
                len(haplotype), 
                haplotype.buffer_size,
                read.qual,
                haplotype.gap_open,
                COMDM.gap_extend,
                COMDM.nucprior,
                read_start_in_hap + firstpos,
                aln1, aln2)

        if ali.contents.score < best_ali_score:
            best_ali_score = ali.contents.score

    return best_ali_score









