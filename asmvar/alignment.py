"""
This module contain the alignment functions for reads to haplotype,
genotyping or realignment process
"""
import os
import sys
import ctypes
import numpy as np
import pysam
import time

import common as com # Same common functions
import datum  as DM  # The global common datum
from read import Read

# Defined some global values
COMDM  = DM.CommonDatum() # DM.CommonDatum().hashmer is defualt to be 7
ALIMER = 15 # Fix by the alignment algorithm in ``align.c``

# Diretory of this python file
dir   = os.path.dirname(os.path.abspath(__file__))
align = ctypes.CDLL(dir + '/align.so') # The alignment module written by C

# Create a class to get the alignment information which record in a 
# c structure data type in `align.fastAlignmentRoutine`
# Thanks for the help of `http://blog.csdn.net/joeblackzqq/article/details/10441017`
class AlignTagPointer(ctypes.Structure):
    _fields_ = [("score", ctypes.c_int), ("pos", ctypes.c_int)]
align.fastAlignmentRoutine.restype = ctypes.POINTER(AlignTagPointer)
###############################################################################

def alignReadToHaplotype(hap_info, bamfile, chrom, sam_idx):
    """ Learn from Platypus.

    This is the most basic and important part of AsmVar. This function decides
    where to anchor the read to the specified haplotype, and calls the banded-
    alignment function -- fastAlignmentRoutine function in 'align.c'. 

    Each bam file represent to one sample. Get the likelihood by 
    reads realignment.

    Args:
        `hap_info`: Element's format is: [start, end, [haplotypes], wvar] 

    Requires pysam
    """
    end_pos = 0
    for h in hap_info:
        if h[1] > end_pos: end_pos = h[1]  

    hapidx2realign_region = {}
    all_indexs = range(len(hap_info))
    pre_start  = None
    bam_reader = pysam.AlignmentFile(bamfile)
    start_loop_idx = 0
    t_n, o_n = 0, 0
    print >> sys.stderr, '[INFO] Now loading the bamfile:', chrom, bamfile
    for al in bam_reader.fetch(chrom):

        alig_start = al.pos + 1 
        alig_end   = al.aend if al.aend else al.pos + al.qlen

        if alig_end > end_pos: break # Don't need to read bam now

        if pre_start and pre_start > alig_start:
            raise ValueError('[ERROR] The bamfile(%s) should be sorted!' % bamfile)
        pre_start = alig_start

        start_loop_idx, overlap_hap_idx = _is_overlap(alig_start, 
                                                      alig_end, 
                                                      start_loop_idx, 
                                                      all_indexs, 
                                                      hap_info)
        is_ovlp = False
        read = Read(al) if overlap_hap_idx else None

        for i in overlap_hap_idx: # It will not be empty if overlap
            
            if i not in hapidx2realign_region:
                regions = _realign_region_in_hap(hap_info[i][2])
                hapidx2realign_region[i] = regions
            # Check whether the region is in the realign region of hap
            loop_idx = 0
            all_idx  = range(len(hapidx2realign_region[i]))
            loop_idx, ovlp_var_idx = _is_overlap(alig_start, 
                                                 alig_end, 
                                                 loop_idx, 
                                                 all_idx,
                                                 hapidx2realign_region[i])
            if ovlp_var_idx:
                # It's Overlap! Now loop all the haplotypes of this 
                # variant block 
                is_ovlp = True
                hap_loglike = []
                for h in hap_info[i][2]: # Loop haplotype 
                    loglikelihood, alis, alie = singleRead2Haplotype(h, 
                                                                     read, 
                                                                     al.pos + 1)
                    # The alignment likelihood of the read map to this hap
                    h.loglikelihood[sam_idx].append(loglikelihood)
                    hap_loglike.append([loglikelihood, alis, alie])

                hap_loglike = np.array(hap_loglike)
                mlhi = hap_loglike[:,0].argmax() # Max likelihood hap idx
                _, alis, alie = hap_loglike[mlhi]
                for j in ovlp_var_idx:
                    # Just assign the read coverage to the haplotype which 
                    # got the max loglikehood
                    if alis < hap_info[i][2][mlhi].variants[j].POS < alie:
                        hap_info[i][2][mlhi].variants[j].cov[sam_idx] += 1

        t_n += 1
        if is_ovlp: o_n += 1
        if t_n % 1000000 == 0:
            print >> sys.stderr, ('[INFO] Loading %d lines and hit POS %d, '
                                  'and %d are happen to overlap. Time: %s' % (
                                   t_n, al.pos + 1, o_n, time.asctime()))

    # Finish looping the whole bamfile of chrom, and we've record all the 
    # likelihood score in haplotype.likelihood of the sample 
    print >> sys.stderr, ('[INFO] Finish Loading %d lines and hit POS %d, '
                          'and %d are happen to overlap. Time: %s\n' % (
                           t_n, al.pos + 1, o_n, time.asctime()))
    bam_reader.close()

def _realign_region_in_hap(haplotypes):
    """
    Find all the variant region in haplotypes

    Args:
        `haplotypes`: It's a list of `Haplotype`
    """
    boundary = 1
    regions  = set()
    for h in haplotypes: # Loop all the haplotypes in this block

        if h.variants:

            for v in h.variants:
                regions.add((max(0, v.POS - boundary), v.POS + boundary))
        else:
            # Small haplotype or it may be a big haplotype but it's a 
            # refernce haplotype and have no variants in it and just 
            # cut the load reads region to fix the max aligment size.
            start = h.hap_start - 1
            end   = min(h.hap_start + boundary, h.hap_end)
    # Cause: Some regions in 'regions' may overlap with others
    # But we don't need to care for this situation!
    return sorted(list(regions))

def _is_overlap(reg_start, reg_end, start_loop_idx, all_idx, regions):
    """
    Chech whether [reg_start, reg_end] overlap with `regions`
    """
    overlap_idx = []
    first_time_ovlp = True
    next_start_idx  = start_loop_idx
    for i in all_idx[start_loop_idx:]:
        
        if reg_start > regions[i][1]: continue
        if reg_end   < regions[i][0]: break
        
        # Overlap now
        if first_time_ovlp:
            next_start_idx  = i
            first_time_ovlp = False
        overlap_idx.append(i)

    return next_start_idx, overlap_idx

def singleRead2Haplotype(haplotype, read, read_align_pos):
    """
    Mapping a single read to a specified haplotype sequence. 

    This is used to find the best anchor position in the haplotype sequence
    for the specified read.

    Args:
        'haplotype': value of ``Haplotype`` class
        'read'     : value of ``Read`` class, which contain read's information 
                     collecting from ``AlignedRead`` of pysam by reading bam
                     files and a hash table for the read sequence

    Return: The best alignement score's likelihood and alignment region
    """
    # hash the read.seqs. And just do it here, and just do it one time!!!
    if not read.seq_hash: # hash the read
        read.seq_hash = com.SeqHashTable(read.seqs, COMDM.hashmer)

    # Firstly, we use seq_hash to find the most possible anchor position
    max_hit = 0
    idx     = {}
    anchor_idx = set() # Record the mapping position in haplotype.sequence
    # Scan read sequence hash table
    for i, id in enumerate(read.seq_hash.hash_pointer): 

        # Ignore the position if hit N!
        if id == COMDM.hitN: continue

        if id in haplotype.seq_hash.hash_table:
            # Index in list of haplotype.seq_hash[id]
            # [BUG] This could be a bug if we hit tandem repeat regions, we may
            # always get the first position at the begin, because we don't use 
            # any strategies to select the most likely start from all the index 
            # in this list(haplotype.seq_hash.hash_table[id])!
            idx[id] = idx.get(id, -1) + 1     # Start from index 0
            if idx[id] >= len(haplotype.seq_hash.hash_table[id]):
                # Rock back! The duplicate hashmer should just allow in 'read'!
                idx[id] -= 1

            pos_idx = haplotype.seq_hash.hash_table[id][idx[id]]
            if pos_idx not in anchor_idx:
                haplotype.map_depth[pos_idx] = 0 # First time -> reset count

            anchor_idx.add(pos_idx)           # Record the mapping pos_idx
            haplotype.map_depth[pos_idx] += 1 # Record the mapping depth
            if haplotype.map_depth[pos_idx] > max_hit:
                max_hit = haplotype.map_depth[pos_idx]

    hap_len_for_align = len(read) + ALIMER
    aln1 = ''.join(['\0' for i in range(2 * len(read) + ALIMER)])
    aln2 = ''.join(['\0' for i in range(2 * len(read) + ALIMER)])
    best_ali_score = 1000000 # A big enough bad score
    best_ali_pos = -1
    firstpos     = 0
    if max_hit > 0:
        # max_hit > 0 means we can find some positions that read could
        # anchor in haplotype by hash searching.

        # Go through all the positions of haplotype.sequence 
        pre_idx = None
        for i in sorted(list(anchor_idx)): 

            if (pre_idx is not None) and (i - pre_idx < 0.5 * len(read)):
                # This will help discount some unnecessary mapping
                continue
            pre_idx = i

            d = haplotype.map_depth[i]
            # Now let's find the best anchor position!
            # 'i' could represent the index of read in haplotype
            is_inside = (i + hap_len_for_align) <= len(haplotype)
            if not is_inside: break

            if d == max_hit and is_inside: 
                # Just start from the most possible position, it's 0-base
                read_start_in_hap = max(0, i - int(round(ALIMER / 2.0)))
                s = read_start_in_hap
                e = s + hap_len_for_align
                # The alignment score is the smaller the better and exactly
                # mapping will be 0, others will larger than 0.
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
                if (COMDM.do_calcu_flank_score and haplotype.buffer_size and
                    ali.contents.score > 0):
                    
                    # This step may cause 'ali.contents.score' to be a negative
                    # value. if negative value be the score, we'll get a 
                    # positive log10 value after times with `COMDM.mot`, and it
                    # means we will get a likelihood > 1.0, could this allow?
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
                    best_ali_pos   = i + 1 # position should be 1-base 

                    # Short-circuit this loop if we find an exact match
                    # (0 is the best score, means exact match)
                    if best_ali_score == 0: 
                        return (0, # log value == 0, means the probability==1.0
                                haplotype.hap_start + best_ali_pos - 1,
                                haplotype.hap_start + best_ali_pos + len(read) - 1)

    # Try original mapping position. If the read is past the end of the 
    # haplotype then don't allow the algorithm to align off the end of 
    # the haplotype. This will only happen in the case of very large deletions.
    read_start_in_hap = min(len(haplotype) - hap_len_for_align, 
                            read_align_pos - haplotype.hap_start)

    if read_start_in_hap != best_ali_pos:
        read_start_in_hap = max(0, i - int(round(ALIMER / 2.0)))
        s = read_start_in_hap 
        e = s + hap_len_for_align 
        # The alignment score is the smaller the better, so that the exactly
        # mapping will be 0, others will larger than 0.
        ali = align.fastAlignmentRoutine(haplotype.sequence[s:e], 
                                         read.seqs, read.qual, 
                                         hap_len_for_align, len(read),
                                         COMDM.gap_extend, 
                                         COMDM.nucprior, 
                                         haplotype.gap_open[s:e],
                                         aln1, aln2)

        if ali.contents.pos != -1:
            firstpos = ali.contents.pos

        if (COMDM.do_calcu_flank_score and haplotype.buffer_size and 
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
            best_ali_pos   = s + 1

    # Finaly, we should convert the score to be a log10 value
    # The probability of read aligne error (convert to be a log10 value)
    log_prob_read_map_error = min(-1e-6, read.mapqual * COMDM.mot)
    # The probability of read aligne correct (still keep the log10 value)
    prob_read_map_right = np.log10(1.0 - 10 ** log_prob_read_map_error)

    log_likelihood_threshold = -100 # A small enough value
    # The max value could just be 0 for all situations if we don't adjust
    # the value with 'do_calcu_flank_score'
    loglk = round(COMDM.mot * best_ali_score + prob_read_map_right, 6)
    return (max(loglk, log_likelihood_threshold),   # It's a log10 value
            haplotype.hap_start + best_ali_pos - 1, # Best align start
            haplotype.hap_start + best_ali_pos + len(read) - 1) # align end









