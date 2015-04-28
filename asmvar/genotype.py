"""
This module contain fucntions and class for genotyping the variats.

This module will connect the `Haplotype` and `alignment` class
"""
import copy
from itertools import combinations as itertools_combinations

import vcf
import alignment as alg # The alignment module
from haplotype import Haplotype

class Diploid(object):
    """
    A class represent a diploid genotype. It store two haplotype
    """
    def __init__(self, haplotype1, haplogtype2):
        """
        Constructor.

        Args:
            `haplotype1`: The data type is Haplotype
            `haplotype2`: The data type is Haplotype
        """
        if hap1.hapregion_hash_id != hap1.hapregion_hash_id:
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

        return a value coresponse to the likelihood.
        """
        # List of likelihood for each aligning read. each value in the array
        # represent a likelihood value of read align to the haplotype
        self.hap1.likelihood = alg.alignReadToHaplotype(self.hap1,
                                                        read_buffer_dict,
                                                        bam_reader)
        self.hap2.likelihood = alg.alignReadToHaplotype(self.hap2,
                                                        read_buffer_dict,
                                                        bam_reader)

def generateAllGenotypes(ref_fa_stream, max_read_len, winvar):
    """
    Generate a list of potentail genotypes.
    """
    gentypes   = []
    haplotypes = _generateAllHaplotypeByVariants(ref_fa_stream, 
                                                 max_read_len, 
                                                 winvar)
    index = range(len(haplotypes))
    for i in index:
        for j in index[i:]:
            # Create Diploid
            gentypes.append(Diploid(haplotypes[i], haplotypes[j]))

    return gentypes

def _generateAllHaplotypeByVariants(ref_fa_stream, max_read_len, winvar):
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
    




