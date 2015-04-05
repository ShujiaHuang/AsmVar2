"""
This module contains various functions and classes for handling variant
information, including utilities for generating combinations of variants, 
haplotypes, and genotypes.
"""
import copy
import logging
import vcf

logger = logging.getLogger('Log')

#####################################################################

class VariantCandidateReader(object):
    """
    A class to read variant information from vcf files, and return
    a batch of stream of variant instances.

    Required PyVCF
    """

    def __init__(self, filenames, options):

        """
        Constructor. Takes the vcf files (bgzip format) and return
        vcf stream for the variant instances.
        """

        self.options     = options
        self.vcf_readers = []

        for filename in filenames:

            if '.gz' not in filename:
                logger.error(
                    '\nSource File %s does not look like a bgzipped VCF '
                    '(i.e. name does not end in .gz)\nYou must supply a '
                    'VCF that has been compressed with bgzip, and indexed '
                    'with tabix' % (filename))
                logger.error(
                    'Compress the VCF using bgzip: '
                    'bgzip %s --> %s.gz' % (filename, filename))
                logger.error(
                    'Remember to use the "tabix -p vcf %s.gz" command to '
                    'index the compressed file' % (filename))
                logger.error(
                    'This should create an index file: %s.gz.tbi\n' % (filename))
                #raise ValueError('\nInput VCF source file %s was not compressed and indexed' % (filename))
            else:
                # How To: how to close the open vcf files' handle?
                self.vcf_readers.append(vcf.Reader(filename = filename))

    def variants(self, chrom = None, start = None, end = None):
        """
        *(chrom = None, start = None, end = None)*
        Generator funtion. Yields variants in order of genomic coordinate.

        fetch variants in a region using 0-based indexing.
        
        The region is specified by :term:`chrom`, *start* and *end*.

        fetch returns an empty string if the region is out of range or
        addresses an unknown *chrom*.

        If *chrom* is given and *start* is None, the sequence from the
        first base is returned. Similarly, if *end* is None, the sequence
        until the last base is returned.

        """
        varlist = []
        # TO DO: All the variant in the same position could make by a net
        for vcf_reader in self.vcf_readers:

            for r in vcf_reader.fetch(chrom, start, end):

                # Continue, if the ALT == '.'
                if r.ALT[0] is None: continue

                # ignore the information that we don't care
                r.INFO    = None
                r.FORMAT  = None
                r.samples = None

                if r.is_snp:
                    varlist.append(r)
                else:
                    # Copy the VCF row, and all the change will just happen 
                    # in this new copy. It would be very useful when we assign
                    # new record to 'varlist'
                    record = copy.copy(r)

                    # TO DO: This trim strategy is greedy, which is not a good 
                    # method. 
                    # e.g: Assume the turth is [A, ATC] 
                    # If we see [AA, ATCA] then will be [A, TCA] instead of
                    # [A, ATC] after triming 
                    for alt in r.ALT:
                        # Non snp variants may leading and/or trailing bases 
                        # trimming.
                        # For indel the first base will always be the same with
                        # REF. That will not be trim!

                        pos, ref = r.POS, r.REF
                        
                        # Trim the leading bases, should keep at lest 1 base 
                        while (alt.sequence[0].upper() == ref[0].upper() and 
                               len(ref) > 1 and len(alt) > 1):
                            alt.sequence = alt.sequence[1:]
                            ref  = ref[1:]
                            pos += 1

                        # Trim the trailing bases, should keep at lest 1 base
                        while (alt.sequence[-1].upper() == ref[-1].upper() and
                               len(ref) > 1 and len(alt) > 1):
                            alt.sequence = alt.sequence[:-1]
                            ref = ref[:-1]

                        if len(ref) > 0 and len(alt) > 0:
                            record.POS = pos
                            record.REF = ref
                            record.ALT = [alt]
                            # After the 'for loop', there may contain some 
                            # duplication positions in 'varlist'
                            varlist.append(record)

        varlist = sorted(list(set(varlist))) # Sorted by reference order
        logger.debug('Found %s variants in region %s in source file' 
                     % (len(varlist), '%s:%s-%s' % (chrom, start, end)))

        # It's a list of '_Record' which type is defined by 'PyVCF'
        return varlist

def get_sequence_context(fa_stream, variant):
    """
    Return the sequence surrounding this variant's position.
    """

    start = max(0, variant.POS - 10)
    return fa_stream.fetch(variant.CHROM, start, variant.POS + 10)

def homoRunForOneVariant(fa_stream, variant):

    """
    Calculate and return the length of the largest homopolymer
    touching this variant. Compute homopolymer lengths on the
    left and right, and return the largest.

    Args:
        fa_stream: A file stream of fastafile, open by pysam.FastaFile
        variant: It's a vcf.model._Record
    """

    left_ref_seq = fa_stream.fetch(variant.CHROM,
                                   max(0, variant.POS - 20),
                                   variant.POS)
    right_ref_seq = fa_stream.fetch(variant.CHROM,
                                    variant.POS,
                                    variant.POS + 20)
    if len(left_ref_seq) == 0 or len(right_ref_seq) == 0:
        return 0

    left_hrun_size  = _calHrunSize(left_ref_seq[::-1].upper())
    right_hrun_size = _calHrunSize(right_ref_seq.upper())

    if left_ref_seq[-1].upper() != right_ref_seq[0].upper():
        return max(left_hrun_size, right_hrun_size)
    else:
        # The left and right side sequence is the same 
        return left_hrun_size + right_hrun_size

def _calHrunSize(sequence):
    """
    Calculate and return the run size of homopolymer

    Args:
        sequence: A string of DNA sequence
    """

    hr = 0
    first_char = sequence[0]
    for c in sequence:
        if c == first_char:
            hr += 1
        else:
            break
    # The hr == 1 means there's not a homopolymer run
    if hr == 1: hr = 0

    return hr
###############################################################################





