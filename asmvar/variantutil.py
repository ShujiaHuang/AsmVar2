"""
This module contains various functions and classes for handling variant
information, including utilities for generating combinations of variants, 
haplotypes, and genotypes.
"""
import copy
import logging
import vcf

import datum as DM  # The global common datum
MODEL = DM.Model()

logger = logging.getLogger('Log')
###############################################################################

class VariantCandidateReader(object):
    """
    A class to read variant information from vcf files, and return
    a batch of stream of variant instances.

    Required PyVCF
    """

    def __init__(self, filenames, options = None):

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
                raise ValueError('\nInput VCF source file %s was not compressed '
                                 'and indexed' % (filename))
            else:
                # How To: how to close the open vcf files' handle?
                self.vcf_readers.append(vcf.Reader(filename = filename))

    def variants(self, chrom = None, start = None, end = None, nosnp = True):
        """
        *(chrom = None, start = None, end = None, nosnp = True)*
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

                if nosnp and r.is_snp:

                    continue
                elif r.is_snp:

                    varlist.append(r)
                else:
                    # Copy the VCF row, and all the change will just happen 
                    # in this new copy. It would be very useful when we assign
                    # new record to 'varlist'
                    record = copy.copy(r)

                    # TO DO: This trim strategy is greedy algorithm, which is 
                    # not the best method. 
                    # e.g: Assume the turth is [A, AATCA] 
                    # If we see [AA, AATCAA] then will be [A, ATCAA] instead of
                    # [A, AATCA] after triming 
                    for alt in r.ALT:
                        # Non snp variants may leading and/or trailing bases 
                        # trimming.
                        # For indel the first base will always be the same with
                        # REF. That will not be trim!

                        pos, ref = r.POS, r.REF
                        # Trim the leading bases, should keep at lest 1 base 
                        while (len(ref) > 1 and len(alt) > 1 and 
                               alt.sequence[:2].upper() == ref[:2].upper()):
                            # At first I think if we just set
                            # `alt.sequence[0].upper() == ref[0].upper()` is
                            # still OK. But after a few seconds, I find that's
                            # wrong! We must always guarrantee the first base
                            # of REF and ALT be the same even after we delete
                            # it. That is why I have to compare the first two
                            # bases insteading of just one!
                            alt.sequence = alt.sequence[1:]
                            ref  = ref[1:]
                            pos += 1

                        # Trim the trailing bases, should keep at lest 1 base
                        while (len(ref) > 1 and len(alt) > 1 and 
                               alt.sequence[-1].upper() == ref[-1].upper()):
                            alt.sequence = alt.sequence[:-1]
                            ref = ref[:-1]

                        if len(ref) > 0 and len(alt) > 0:
                            record.POS = pos
                            record.REF = ref
                            record.ALT = [alt]
                            # After the 'for loop', there may contain some 
                            # duplication positions in 'varlist'
                            varlist.append(record)

        varlist = self._dedup(varlist)
        logger.debug('Found %s variants in region %s in source file' 
                     % (len(varlist), '%s:%s-%s' % (chrom, start, end)))

        # It's a list of '_Record' which type is defined by 'PyVCF'
        return varlist

    def _dedup(self, varlist):
        """
        Delete conflict variants at the same position.  

        And keep the small one
        """
        if len(varlist) == 0: 
            return [] 

        vdict = {}
        # Firstly, we should the duplicate positions
        for i, v in enumerate(varlist):
            k = v.CHROM + ':' + str(v.POS)
            vdict[k] = vdict.get(k, []) + [i]

        variants = []
        for k, idx in vdict.items():

            if len(idx) == 1:
                variants.append(varlist[idx[0]])
            else:
                """
                delete conflict here
                """
                prevar = varlist[idx[0]]
                for i in idx[1:]:

                    if prevar.is_snp and (not varlist[i].is_snp):
                        # Ignore SNP,if we find other variant type  
                        prevar = varlist[i]
                    elif prevar.REF == varlist[i].REF:
                        # Merge variants if they have the same REF, SNP
                        # will be merge here, too
                        prevar.ALT += varlist[i].ALT
                    elif varlist[i].is_snp:
                        # Ignore SNP,if we find other variant type
                        continue
                    else:
                        # Keep the smaller one
                        pre_max_vlen  = max([abs(len(prevar.REF) - len(v)) 
                                             for v in prevar.ALT])
                        this_max_vlen = max([abs(len(varlist[i].REF) - len(v)) 
                                             for v in varlist[i].ALT])
                        if this_max_vlen < pre_max_vlen:
                            prevar = varlist[i]
                        
                prevar.ALT = list(set(prevar.ALT)) # Unique the sequence
                variants.append(prevar)

        return sorted(list(set(variants))) # Sorted by reference pos order 

def calPrior(fa_stream, variant):
    """
    calculate and return the prior probability for this variant.

    Use the module of Pindel

    Args:
        'variant': It's variant record by module: 'vcf.Reader'
    """
    prior = 0.0

    if variant.is_snp:
        # SNP
        prior = 1e-3 / 3

    elif len(variant.REF) == len(variant.ALT[0]):
        # Substitution
        diff_num = len(
            [(x, y) for (x, y) in zip(variant.REF, variant.ALT[0]) if x != y])
        prior = 5e-5 * (0.1 ** (diff_num - 1)) * (1.0 - 0.1)

    elif len(variant.REF) == 1:
        # Insertion. Use the most easy model this moment
        # And I'll update it later as 'Platypus' in 'variant.calculatePrior'
        prior = 1e-4 * 0.25 ** (len(variant.ALT[0]) - 1)
    elif len(variant.ALT[0]) == 1: 
        # Deletion. Use the most easy model this moment
        # And I'll update it later as 'Platypus' in 'variant.calculatePrior'
        prior = 1e-4 * 0.6 ** (len(variant.REF) - 1)
    else:
        # Replacement 
        prior = 5e-6

    return max(1e-10, prior)

#def _indelPrior(fa_stream, variant, indel_size):
    #"""
    #Calculate indel prior, based on sequence context.
    #"""
    # MODEL
    # context_size = 100
    # sequence = get_sequence_context(fa_stream, variant)

def get_sequence_context(fa_stream, variant, size = 10):
    """
    Return the sequence surrounding this variant's position.
    """

    start = max(0, variant.POS - size)
    return fa_stream.fetch(variant.CHROM, start, variant.POS + size)

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





