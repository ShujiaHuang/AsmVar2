"""
This module contain fucntions and class for genotyping the variats.

This module will connect the `Haplotype` and `alignment` class
"""
import copy
from itertools import combinations as itertools_combinations
import numpy as np

import vcf
import alignment as alg # The alignment module
from haplotype import Haplotype

class Diploid(object):
    """
    A class represent a diploid genotype. It store two haplotype
    """
    def __init__(self, haplotype1, haplotype2):
        """
        Constructor.

        Args:
            `haplotype1`: The data type is Haplotype
            `haplotype2`: The data type is Haplotype
        """
        if haplotype1.hapregion_hash_id != haplotype2.hapregion_hash_id:
            raise ValueError('[ERROR] The two haplotypes did not have the same '
                             'haplotype region, could not build up a diploid.')

        self.hap1 = haplotype1 # Do I have use deepcopy here?
        self.hap2 = haplotype2 # Do I have use deepcopy here?

    def calLikelihood(self, read_buffer_dict, bam_reader):
        """
        Calculate the genotype likelihood for this single bam_reader's sample.

        Args:
            `read_buffer_dict`: A hash to reads. It's a contianer of reads, 
                                use it to prevent recalculting the hash-sequence 
                                for the same reads. This could save the running 
                                time.
            `bam_readers`: A single bamfile reader opened by `pysam.AlignmentFile`

        return a log10 likelihood.
        """
        # List of likelihood for each aligning read. each value in the array
        # represent a likelihood value of read align to the haplotype
        # Remember these likelihood will be changed follew different bamfile
        self.hap1.likelihood = alg.alignReadToHaplotype(self.hap1,
                                                        read_buffer_dict,
                                                        bam_reader)
        self.hap2.likelihood = alg.alignReadToHaplotype(self.hap2,
                                                        read_buffer_dict,
                                                        bam_reader)
        lksize1 = len(self.hap1.likelihood)
        lksize2 = len(self.hap2.likelihood)
        if lksize1 != lksize2:
            raise ValueError('[ERROR] The two haplotype in this diploid must'
                             'have the same number alignment of reads. But '
                             'the number is %s and %s .' % (lksize1, lksize2))

        likelihood = None
        if lksize1: # Not empty 
            likelihood = 0.0

        # Calculate the log10 likelihood for this diploid region
        for i in range(lksize1):

            log10lk1 = self.hap1.likelihood[i] 
            log10lk2 = self.hap2.likelihood[i] 

            # Almost good to 1000 times. Just take the highest and forget 
            # the small one
            if abs(log10lk1 - log10lk2) >= 3.0:
                likelihood += (np.log10(0.5) + max(log10lk1, log10lk2))

            # The likelihoods of the two haplotypes are almost the same. This
            # could just happens when either they both bad fits or both are 
            # perfect match. Any way the read lies in a position that it 
            # cannot be used to distingush the two haplotypes.
            elif abs(log10lk1 - log10lk2) <= 1e-3:
                likelihood += log10lk1

            # Calculate as the normal way: Combine the two likelihood
            else:
                prob = 0.5 * (np.power(10, log10lk1) + np.power(10, log10lk2))
                likelihood += np.log10(prob) # now it's a log10 value

        # Generally, for our situation the likelihood should always < 0, but 
        # it may > 0, once we adjust with flank region of haplotype in the 
        # `alignReadToHaplotype` process.
        read_map_count = lksize1 # Count of mapping reads of this sample

        return likelihood, read_map_count


def generateAllGenotypes(haplotypes):
    """
    Generate a list of potentail genotypes. 
    """
    gentypes   = []
    index = range(len(haplotypes))
    for i in index:
        for j in index[i:]:
            # Create Diploid
            # Each element is [diploid, hap_hash_id1, hap_hash_id2]
            # hash_id is the identify code of each haplotype!
            gentypes.append([Diploid(haplotypes[i], haplotypes[j]), 
                             hash(haplotypes[i]), hash(haplotypes[j])])

    return gentypes

def generateAllHaplotypeByVariants(ref_fa_stream, max_read_len, winvar):
    """
    Generate all the potential haplotype in the window of `winwar`.
    And return a list of haplotype, corresponing to all the possible
    combination of variants in this specific window

    Args:
        `ref_fa_stream`: Input fasta stream of reference
        `max_read_len`: The max length of reads
        `winvar`: It's dict(chrom=chrom, start=start, end=end, variant=var)
    """

    # This is the reference haplotype attempt that no variants in this window
    refhap = Haplotype(ref_fa_stream, winvar['chrom'], winvar['start'], 
                       winvar['end'], max_read_len)
    haplist = []
    num_var = len(winvar['variant'])

    # Generate all the combination of haplotypes
    # All the combinantion of haplotype may be too much! 
    for n in range(num_var):
        n += 1
        for varlist in itertools_combinations(winvar['variant'], n):

            # Must be deep copy, prevent to changing the raw `winvar`
            varlist   = copy.deepcopy(varlist)
            need_iter = True
            while need_iter:
                """
                Iter all the possible variant combination, and each
                variant combination could be join to be one haplotype
                """
                need_iter = False
                var = []
                for i, v in enumerate(varlist):

                    tmp = copy.deepcopy(v) # Still have to deep copy
                    if len(v.ALT) > 1:
                        # v.ALT array size will be smaller here
                        need_iter = True #  
                        tmp.ALT   = [v.ALT.pop()]
                    var.append(copy.deepcopy(tmp))

                hap = Haplotype(ref_fa_stream, winvar['chrom'], winvar['start'],
                                winvar['end'], max_read_len, var)
                haplist.append(hap)

    return haplist




