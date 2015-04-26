"""
This module caontain fucntions and class for genotyping the variats.
"""
import copy
from itertools import combinations as itertools_combinations

import vcf
from haplotype import Haplotype

class Diploid(object):
    """
    A class represent a diploid genotype. It store two haplotype
    """
    def __init__(self, Haplotype hap1, Haplotype hap2):
        """
        Constructor.
        """
        self.hap1 = hap1 # Do I have use deepcopy here?
        self.hap2 = hap2 # Do I have use deepcopy here?

    del likelihood(self, bam_readers):
        """
        Calculate the likelihood for this diploid.

        Args:
            `bam_readers`: A list of bamfiles' input stream
        """

def generateAllGenotypes(ref_fa_stream, max_read_len, winvar):
    """
    Generate a list of potentail genotypes.
    """
    gentypes   = []
    haplotypes = _generateAllHaplotypeByVariants(ref_fa_stream, max_read_len, 
                                                 winvar)
    index = range(len(haplotypes))
    for i in index:
        for j in index[i:]:
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
    




