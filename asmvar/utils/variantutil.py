"""
This module contains various functions and classes for handling variant
information, including utilities for generating combinations of variants, 
haplotypes, and genotypes.
"""

import copy
import logging
import vcf

logger = logging.getLogger('Log')
###############################################################################

class VariantCandidateReader(object):
    """
    A class to read variant information from vcf files, and return
    a batch of stream of variant instances.

    Required PyVCF
    """

    def __init__(self, filename, options = None):

        """
        Constructor. Takes the vcf files (bgzip format) and return
        vcf stream for the variant instances.
        """

        self.options = options
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
            self.vcf_reader = vcf.Reader(filename = filename)


class Variant(object):
    """
    A class to record VCF variants
    """
    def __init__(self, r):
        """
        Generator funtion. Yields variants in order of genomic coordinate.

        fetch variants in a region using 0-based indexing.
        
        The region is specified by :term:`chrom`, *start* and *end*.

        fetch returns an empty string if the region is out of range or
        addresses an unknown *chrom*.

        If *chrom* is given and *start* is None, the sequence from the
        first base is returned. Similarly, if *end* is None, the sequence
        until the last base is returned.

        Args:
            `done_load_var`: A set for recording all the have loaded variants
                             The reason that why I have to use this parameters,
                             it's we may load some variants multiple times if 
                             the variants' region is so big, and the start-end
                             region can not hold them all! The parameter is 
                             'set()' type in python, and I use it to ignore the
                             loaded variants!
        """
        # ignore the information that we don't care
        r.INFO    = None
        r.FORMAT  = None
        r.samples = None
        r.cov     = [] # Add a new value to record the coverage
        self.record = r # vcf.model._Record type
        self.hrun   = None # Size of homo run around the variant
        self.nratio = None # N ratio around the variant

    def parse(self):
        """
        Parse variant and return variant list
        """
        # Copy the VCF row, and all the change will just happen 
        # in this new copy. It would be very useful when we assign
        # new record to 'varset'
        record = copy.deepcopy(self.record)
        varset = set()

        # Do not trim ALT.sequence if ALT is more than one
        #if len(r.ALT) > 1:
        #    varset.add(record)
        #    continue

        # TO DO: This trim strategy is greedy algorithm, which is 
        # not the best method. 
        # e.g: Assume the turth is [A, AATCA] 
        # If we see [AA, AATCAA] then will be [A, ATCAA] instead of
        # [A, AATCA] after triming 
        for alt in self.record.ALT:
            # Non snp variants may leading and/or trailing bases 
            # trimming.
            # For indel the first base will always be the same with
            # REF. That will not be trim!

            pos, ref = self.record.POS, self.record.REF
            rh = 0 # header index of Ref-seq 0-base
            ah = 0 # header index of Alt-seq 0-base
            rt = len(ref) # tail index of Ref-seq 
            at = len(alt) # tail index of Alt-seq
            # Trim the leading bases, should keep at lest 1 base 
            while ((rt - rh > 1) and (at - ah > 1) and 
                (alt.sequence[ah:ah+2].upper() == ref[rh:rh+2].upper())):
                # At first I think if we just set
                # `alt.sequence[0].upper() == ref[0].upper()` is
                # still OK. But after a few seconds, I find that's
                # wrong! We must always guarrantee the first base
                # of REF and ALT be the same even after we delete
                # it. That is why I have to compare the first two
                # bases insteading of just one!
                rh  += 1
                ah  += 1
                pos += 1

            # Trim the trailing bases, should keep at lest 1 base
            while (rt - rh > 1 and at - ah > 1 and 
                   alt.sequence[at-1].upper() == ref[rt-1].upper()):
                rt -= 1
                at -= 1

            if rt - rh > 0 and at - ah > 0:
                record.POS = pos
                record.REF = ref[rh:rt]
                alt_seq    = alt.sequence[ah:at]
                record.ALT = [vcf.model._Substitution(alt_seq)]
                # After this 'for loop', there may contain some 
                # position duplication and overlap positions in 
                # 'varset'
                varset.add(record)

        # Sorted by reference pos order  
        varlist = self._dedup(sorted(list(varset))) 

        # It's a list of '_Record' which type is defined by 'PyVCF'
        return varlist

    def _dedup(self, varlist):
        """
        Delete conflict variants at the same position.  

        And keep the small one
        """
        if len(varlist) == 0: 
            return [] 

        # Firstly, we should detect the duplicate positions
        vdict, max_vlen = {}, {}
        for i, v in enumerate(varlist):
            maxlen = max([abs(len(v.REF) - len(a)) for a in v.ALT])
            k      = v.CHROM + ':' + str(v.POS)
            vdict.setdefault(k, []).append(i)
            max_vlen.setdefault(k, []).append(maxlen)

        variants = set()
        for k, idx in vdict.items():

            if len(idx) == 1:
                variants.add(varlist[idx[0]])

            else:
                """
                delete conflict here
                """
                prevar, pre_max_vlen = varlist[idx[0]], max_vlen[k][0]
                for j, i in enumerate(idx[1:]):
                    j += 1 # Because idx array start from the second one!

                    # Max variant size in varlist[i]
                    if prevar.is_snp and (not varlist[i].is_snp):
                        # Ignore SNP,if we find other variant type  
                        prevar = varlist[i]

                    elif prevar.REF == varlist[i].REF:
                        # Merge variants if they have the same REF, SNP
                        # will be merge here, too
                        for av in varlist[i].ALT:
                            # Avoid deplicate
                            if av not in prevar.ALT:
                                prevar.ALT.append(av)

                    elif max_vlen[k][j] == 0 or varlist[i].is_snp:
                        # Ignore SNP,if we find other variant type
                        continue
                    else:
                        # Keep the smaller one and reset the variant's size
                        if max_vlen[k][j] < pre_max_vlen:
                            prevar       = varlist[i]
                            pre_max_vlen = max_vlen[k][j] # Must re-assign here

                    # Re-assign the max variant size in 'prevar'
                    if pre_max_vlen < max_vlen[k][j]:
                        pre_max_vlen = max_vlen[k][j]

                variants.add(prevar)

        return list(variants)


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
            [(x, y) for (x, y) in zip(variant.REF, str(variant.ALT[0])) if x != y])
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

#def _indelPrior(chr_fa_seq, variant, indel_size):
    #"""
    #Calculate indel prior, based on sequence context.
    #"""
    # MODEL
    # context_size = 100
    # sequence = get_sequence_context(chr_fa_seq, variant)

def get_sequence_context(chr_fa_seq, variant, size = 10):
    """
    Return the sequence surrounding this variant's position.
    """

    start = max(0, variant.POS - size)
    return chr_fa_seq[start:variant.POS + size]

def nRatioForOneVariant(chr_fa_seq, variant, size = 10):
    """
    Calculate the N ratio surrounding this variant's position 
    """
    left_ref_seq  = chr_fa_seq[max(0, variant.POS - size):max(variant.POS - 1, 0)]
    right_ref_seq = chr_fa_seq[variant.POS:variant.POS + size]  
    nr = 0.0
    for v in variant.ALT:
        seq = left_ref_seq + str(v) + right_ref_seq
        n = seq.upper().count('N') / float(len(seq))
        if n > nr:  # Just get the biggest n-ratio
            nr = round(n, 3)

    return nr

def homoRunForOneVariant(chr_fa_seq, variant):
    """
    Calculate and return the length of the largest homopolymer
    touching this variant. Compute homopolymer lengths on the
    left and right, and return the largest.

    Args:
        chr_fa_seq: A fasta sequence of chromosome in `varianti.CHROM`
        variant: It's a vcf.model._Record
    """

    left_ref_seq  = chr_fa_seq[max(0, variant.POS - 20):
                               max(variant.POS - 1, 0)].upper()
    right_ref_seq = chr_fa_seq[variant.POS:variant.POS + 20].upper()
    if len(left_ref_seq) == 0 or len(right_ref_seq) == 0:
        return 0

    left_hrun_size  = _calHrunSize(left_ref_seq[::-1])
    right_hrun_size = _calHrunSize(right_ref_seq)

    if left_ref_seq[-1] != right_ref_seq[0]:
        return max(left_hrun_size, right_hrun_size)
    else:
        # The left and right side sequence is the same 
        if variant.REF[0].upper() == right_ref_seq[0]:
            right_hrun_size += 1
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





