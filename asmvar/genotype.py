"""
This module contain fucntions and class for genotyping the variats.

This module will connect the `Haplotype` and `alignment` class
"""
import sys
import copy
from itertools import combinations as itertools_combinations
import numpy as np

import time

import vcf
import alignment as alg # The alignment module
from haplotype import Haplotype

import datum as DM
COMDM = DM.CommonDatum()

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

        self.hap1 = haplotype1 # Do not modify any thing in `self.hap1`
        self.hap2 = haplotype2 # Do not modify any thing in `self.hap2`

        # The size of this likelihood is depandent by the sample number
        # and it's the same with the row size of self.hap1.loglikelihood or 
        # self.hap2.loglikelihood. And the likelihood in `self.likelihood` is
        # not the log value any more, it's a probability now.
        # pop likelihood means the likelihood of population
        self.pop_loglikelihood = self._set_genotype_likelihood() 

    def _set_genotype_likelihood(self):
        """
        The genotype likelihood must be 1D array
        Calculate the genotype likelihood for this single bam_reader's sample.

        return a log10 likelihood.
        """
        # The column is a list of likelihood for each aligning read. 
        # each row represent to a sample.
        s_size1 = len(self.hap1.loglikelihood)
        s_size2 = len(self.hap2.loglikelihood)
        if s_size1 != s_size2:
            raise ValueError('[ERROR] The two haplotype in this diploid must'
                             'have the same number alignment of sample.  But'
                             'the number is %s and %s.' % (s_size1, s_size2))

        loglikelihoods = []
        for i in range(s_size1):
            # loop sample
            loglikelihoods.append(self.__calIndividualLikelihood(
                self.hap1.loglikelihood[i], self.hap2.loglikelihood[i]))
        return loglikelihoods

    def __calIndividualLikelihood(self, loglikelihood1, loglikelihood2):

        lksize1 = len(loglikelihood1)
        lksize2 = len(loglikelihood2)
        if lksize1 != lksize2:
            raise ValueError('[ERROR]The two haplotype in this diploid must '
                             'have the same number alignment of reads.  But '
                             'the number is %s and %s.' % (lksize1, lksize2))

        # Log value
        likelihood = 0.0 if lksize1 else None # May be no read coverage
        for i in range(lksize1):

            log10lk1 = loglikelihood1[i] 
            log10lk2 = loglikelihood2[i] 

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
                prob = 0.5 * (10 ** log10lk1 + 10 ** log10lk2)
                likelihood += np.log10(prob) # now it's a log10 value

        # Generally, for our situation the likelihood should always < 0, but 
        # it may > 0, once we adjust with flank region of haplotype in the 
        # `alignReadToHaplotype` process.
        best_loglikelihood = -1e20
        return max(best_loglikelihood, likelihood) # log value

def generateAllGenotypes(haplotypes):
    """
    Generate a list of potentail genotypes. 
    """
    gentypes = []
    index = range(len(haplotypes))
    for i in index:
        for j in index[i:]:
            # Create Diploid
            gentypes.append(Diploid(haplotypes[i], haplotypes[j]))

    return gentypes

def generateAllHaplotypeByVariants(chr_fa_seq, max_read_len, winvar):
    """
    Generate all the potential haplotype in the window of `winwar`.
    And return a list of haplotype, corresponing to all the possible
    combination of variants in this specific window

    Args:
        `chr_fa_seq`: The fasta sequence of chromsome in `winvar`
        `max_read_len`: The max length of reads
        `winvar`: It's dict(chrom=chrom, start=start, end=end, variant=var)
    """
    # This is the reference haplotype attempt that no variants in this window
    #refhap = Haplotype(chr_fa_seq, winvar['chrom'], winvar['start'], 
    #                   winvar['end'], max_read_len)
    haplist = [] # Initial with empty list
    num_var = len(winvar['variant'])
    vindex  = range(num_var) # An array contain all the index of variants
    var_overlap_index = _get_overlap_var_index(winvar['variant'])

    # Generate all the combination of haplotypes
    # All the combinantion of haplotype may be too much! 
    for idxs in itertools_combinations(vindex, num_var):
        # Because we have REF in ALT array, swe do not have to iter from 1 but 
        # from 'num_var' directly

        new_idxs = _recreate_varaint_idxlist_by_overlop(idxs, 
                                                        winvar['variant'],
                                                        var_overlap_index)
        for idx in new_idxs:

            # Must be deep copy, prevent to changing the raw `winvar`
            varlist = [copy.deepcopy(winvar['variant'][j]) for j in idx]
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

                hap = Haplotype(chr_fa_seq, winvar['chrom'],
                                winvar['start'], winvar['end'], 
                                max_read_len, var)
                haplist.append(hap)

    return haplist

def _get_overlap_var_index(varants):
    """
    Self overlap. Find all the overlap relationship with other element
    in 'variants' and return a hash which record the overlap index for
    all the array's elements
    """
    var_overlap_index = {}
    all_indexs = range(len(varants))
    loop_i     = 1
    for idx in all_indexs:

        start = varants[idx].POS
        end   = varants[idx].POS + len(varants[idx].REF) - 1
        flag  = True
        var_overlap_index[idx] = [idx] # Initial with the own index

        for i in all_indexs[loop_i:]:
            if idx == i: continue
            if start > varants[i].POS + len(varants[i].REF) - 1: continue
            if end < varants[i].POS: break

            # Overlap now!
            if flag:
                loop_i = i
                flag   = False

            var_overlap_index[idx].append(i)

    return var_overlap_index

def _recreate_varaint_idxlist_by_overlop(index_list, 
                                         variants,
                                         index_overlap_hash): 
    """
    Separating the overlap variants' index into different haplotype.
    
    return the new index list and it's a 2-D array.
    ####
    Args:
        `index_list`: A list of array's indexs
        `index_overlap_hash`: A hash contain the overlap relationship
    """
    pre_index   = index_list[0]
    index_nodes = [[pre_index]] # The first index
    for i in index_list[1:]:
        # Getting Overlap and creat index nodes
        if i in index_overlap_hash[pre_index]:
            # Collecte all the overlap variant's index into one node 
            index_nodes[-1].append(i)
        else:
            index_nodes.append([i]) # New node

        if (variants[pre_index].POS + len(variants[pre_index].REF) - 1 < 
            variants[i].POS + len(variants[i].REF) - 1):
            pre_index = i

    new_index_list = []
    need_iter      = True
    while need_iter:

        need_iter = False
        tmp_index_list = []
        for index in index_nodes:
            # get index by looping all the index nodes.
            i = index[0] # default
            if len(index) > 1:
                need_iter = True
                i = index.pop()

            tmp_index_list.append(i)

        new_index_list.append(tmp_index_list)

    return new_index_list

