"""
This module contain class for Diploid.
"""
import sys
import numpy as np

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

