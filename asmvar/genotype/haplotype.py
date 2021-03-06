"""
This module contain classes and functions for general haplotype.
"""
import sys
import copy
import logging

sys.path.append('..')
from utils import common as com # Same common functions
from utils import datum  as DM  # The global common datum
from utils import variantutil as vutil

COMDM  = DM.CommonDatum() # DM.CommonDatum().hashmer is defualt to be 7
logger = logging.getLogger('Log')

class Haplotype(object):
    """
    Class to encapsulate a single haplotype from a list of variants and 
    reference region 'chrId, start, end' fasta file's stream.

    It'll just return the reference haplotype if `variants` is None

    Args:
        chrom: reference chromosome id
        start_pos: reference start position (1-based coordinate systems)
        end_pos: reference end position
        chr_fa_seq: The fasta sequence of `chrom`.
        max_read_length: The max size of reads in fq file, use for extending 
                         the haplotype region. [100]
        variants: A list of variants. The data type defined in `PyVCF` module
                  and been called by `VariantCandidateReader` in the module 
                  `variantutil`. [None]
    """
    def __init__(self, chr_fa_seq, chrom, start_pos, end_pos, 
                 max_read_length, variants = []):
        """
        Initial the Haplotype
        """
        chromlen         = len(chr_fa_seq)
        self.chrom       = chrom
        self.buffer_size = min(150, 1.5 * max_read_length) # a reasonable size
        self.start_pos   = max(1, start_pos)             # Window start 
        self.end_pos     = min(end_pos, chromlen)        # window end
        self.hap_start   = max(1, self.start_pos - self.buffer_size)
        self.hap_end     = min(self.end_pos + self.buffer_size, chromlen)
        self.variants    = copy.deepcopy(variants) # Must use deepcopy!!
        self.sequence    = None # The haplotype sequence
        self.hash_id     = None # A hash id use for identified this haplotype

        # This is a haplotype region identify. and could be used for guarantee 
        # the two haplotypes have the same region when we construct diploid.
        self.hapregion_hash_id = hash((self.chrom, self.start_pos, self.end_pos))

        if len(variants) == 0:
            self.sequence = chr_fa_seq[self.hap_start - 1:self.hap_end]
        else:
            # Left side of variants, donot include self.start_pos base
            start1, end1 = self.hap_start, max(self.start_pos - 1, 1)
            # Right side of variants, donot include self.end_pos base
            start2, end2 = self.end_pos, self.hap_end
            
            # CAUTION: 'start1' and 'start2' treat as 0-base
            leftseq  = chr_fa_seq[start1:end1]
            rightseq = chr_fa_seq[start2:end2]
            self.sequence = (leftseq + 
                             self._getMutatedSequence(chr_fa_seq) + 
                             rightseq)

        # Record the mapping depth for each position of haplotype
        self.map_depth = [0] * len(self.sequence)
        # Record the score penalize of gap open. Use it for read realign process
        # The size of ``gap_open`` should be the same as self.sequence, but
        # now I'll assign it to be None, and I'll assign a string for it when
        # we need it and the ASCII of each charter of the is a penalize value
		# for gap open
        self.gap_open = com.set_gap_open_penalty(self.sequence, 
                                                 COMDM.homopol_penalty)
        # Encode a hash table for self.sequence used for mapping when we needed
        self.seq_hash = com.SeqHashTable(self.sequence, COMDM.hashmer)
        # loglikelihood calculated for each read during re-aligning process
        self.loglikelihood = []

    def __len__(self):
        return len(self.sequence) # The sequence's size

    def __hash__(self):
        """
        This function allows haplotypes to be hashed, and so stored in a set 
        or dictionary. The supporting reads are not included in the hashing, 
        as we want two haplotypes to give the same hash id if they have the 
        same positions and sequences.
        """
        if self.hash_id is None:
            self.hash_id = hash((self.chrom, self.start_pos, self.end_pos, 
                                 self.sequence))

        return self.hash_id

    def homoRunLength(self, chr_fa_seq):
        """
        I don't think we need this function here!
        """
        return [vutil.homoRunForOneVariant(chr_fa_seq, v) for v in self.variants]

    def _getMutatedSequence(self, fa_sequence):
        """
        Return the sequence mutated with all the variants being considered
        in this haplotype

        Just to remind ourselves: SNPs are reported at the index of the changed
        base (0-indexed internally, and 1-indexed in VCF). Insertions are reported
        such that the index is that of the last reference base before the insertion.
        Deletions are reported such that the index is that of the first deleted base.
        """
        if len(self.variants) == 0:

            region = '%s:%d-%d' % (self.chrom, self.start_pos, self.end_pos)
            logger.error('Region(%s) does not contain variants.' % region)
            raise ValueError('Region(%s) does not contain variants. Do not '
                             'have to call _getMutatedSequence' % region)
        seq = []
        current_pos = self.start_pos

        for v in self.variants:

            if len(v.ALT) != 1:
                # For each haplotype the ALT should just be one.
                # It should be divide into > 1 haplotypes before we call this 
                # function to create the hap-mutate-seq if len(v.ALT) > 1, so
                # duplication positions should be NEVER happen!
                raise ValueError('The ALT sequence should just be one '
                                 'for each single haplotype!')

            if current_pos > v.POS: 
                # Deplication position should never happen!
                raise ValueError('Start pos(%d) > end pos(%d). The variant '
                                 'may not sorted or duplicate or overlap with '
                                 'the preview variant' % (current_pos, v.POS))

            # current_pos is always point to the last breakpoint
            if current_pos == v.POS and v.ALT[0] is not None:
                seq.append(v.ALT[0].sequence)
            else: # current_pos < v.POS
                # Get sequence up to one base before the variant. And DO NOT
                # contain the refernce base on v.POS, which one is already
                # the first base of mutate-seq on v.POS
                ref = fa_sequence[current_pos:v.POS - 1]
                # Mutate-seq v.POS
                if v.ALT[0] is not None:
                    seq.append(ref + v.ALT[0].sequence)
                else:
                    seq.append(ref)
            # Move up to the next variant
            current_pos = v.POS + len(v.REF) - 1    # 0-base!

        # Now we should deal with the last postion
        if current_pos > self.end_pos:
            raise ValueError('Start pos(%d) > end pos(%d). The variant may not '
                             'sorted or duplicate or overlap with the preview '
                             'variant...' % (current_pos, self.end_pos))
            
        if current_pos < self.end_pos:
            seq.append(fa_sequence[current_pos:self.end_pos])
        # Join them to be a single string
        return ''.join(seq)


