"""
This module contain classes and functions for general haplotype
"""
import copy
import logging

import variantutil as vutil

logger = logging.getLogger('Log')

class Haplotype(object):
    """
    Class to encapsulate a single haplotype from a list of variants and 
    reference region 'chrId, start, end' fasta file's stream.

    It'll just return the reference haplotype if ``variants`` is None

    Args:
        chrom: reference chromosome id
        start_pos: reference start position (1-based coordinate systems)
        end_pos: reference end position
        ref_stream_fa: A file stream of reference fasta file.
        max_read_length: The max size of reads in fq file, use for extending 
                         the haplotype region. [100]
        variants: A list of variants. The data type defined in ``vcf`` module
                  and been called by ``VariantCandidateReader`` in the module 
                  ``variantutil``. [None]
    """
    def __init__(self, ref_stream_fa, chrom, start_pos, end_pos, 
                 max_read_length = 100, variants = None):

        chromlen         = ref_stream_fa.get_reference_length(chrom)

        self.chrom       = chrom
        self.buffer_size = min(500, 2 * max_read_length) # a reasonable size
        self.start_pos   = max(1, start_pos)             # Window start 
        self.end_pos     = min(end_pos, chromlen)        # window end
        self.hap_start   = max(1, self.start_pos - self.buffer_size)
        self.hap_end     = min(self.end_pos + self.buffer_size, chromlen)
        self.fa_stream   = ref_stream_fa
        self.variants    = copy.copy(variants) # Must use copy!!
        self.sequence    = None # The haplotype sequence

        exdstart = max(1, self.start_pos - self.buffer_size) # extend start
        exdend   = min(self.end_pos + self.buffer_size, chromlen) # extend end
        if variants is None or len(variants) == 0:

            self.sequence = ref_stream_fa.fetch(self.chrom,
                                                self.hap_start - 1,
                                                self.hap_end)
        else:
            # Left side of variants, donot include self.start_pos base
            start1, end1 = self.hap_start, max(self.start_pos - 1, 1)
            # Right side of variants, donot include self.end_pos base
            start2, end2 = self.end_pos, self.hap_end
            
            # CAUTION: 'start1' and 'start2' treat as 0-base
            leftseq       = ref_stream_fa.fetch(self.chrom, start1, end1)
            rightseq      = ref_stream_fa.fetch(self.chrom, start2, end2)
            self.sequence = leftseq + self._getMutatedSequence() + rightseq

    def homoRunLength(self):
        """
        I don't think we need this function here!
        """
        return [vutil.homoRunForOneVariant(self.fa_stream, v) for v in self.variants]

    def _getMutatedSequence(self):
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
				# deplication positions should be NEVER happen!
                raise ValueError('The ALT sequence should just be one '
                                 'for each single haplotype!')

            if current_pos > v.POS: 
                # Deplication position should never happen!
                raise ValueError('Start pos(%d) > end pos(%d). The variant '
                                 'may not sorted or duplicate or overlap with '
                                 'the reference region' % (current_pos, v.POS))

            # current_pos is always point to last breakpoint
            if current_pos == v.POS:
                seq.append(v.ALT[0].sequence)
            else: # current_pos < v.POS
                # Get sequence up to one base before the variant. And DO NOT
                # contain the refernce base on v.POS, which one is already
                # the first base of mutate-seq on v.POS
                ref = self.fa_stream.fetch(self.chrom, current_pos, v.POS - 1)
                # Mutate-seq v.POS
                seq.append(ref + v.ALT[0].sequence)
            # Move up to the next variant
            current_pos = v.POS + len(v.REF) - 1    # 0-base!

        # Now we should deal with the last postion
        if current_pos > self.end_pos:
            raise ValueError('Start pos(%d) > end pos(%d). The variant may not '
                             'sorted or duplicate or overlap with the reference'
                             'region.' % (current_pos, self.end_pos))
            
        if current_pos < self.end_pos:
            seq.append(self.fa_stream.fetch(self.chrom, current_pos, self.end_pos))
        # Join into a single string
        return ''.join(seq)
        
