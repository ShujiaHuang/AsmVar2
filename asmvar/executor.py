"""
This module will contain all the executor steps of variant calling,
genotying, variants' recalibration et.el

We have many important modules in AsmVar while this module is the lord
to rule them all, in a word, it's "The Ring". 
`Asmvar.py` is "Sauron", and this module could just be called by it.
"""
import pysam
import vcf
import numpy as np

import variantutil as vutil
import genotype as gnt

import datum as DM  # The global common datum 
COMDM = DM.CommonDatum()

class VariantCaller(object):
    """
    A class for variants' calling.

    Input required:
        (1) Long reads or genome VS genome alignment files, sam/bam format
            and must be sorted
        (2) Target genome(Reference) fasta, has been index (.fai)
        (3) Query genome fasta, has been index (.fai)
    """
    def __init__(self, alignfile, targetfile, queryfile, options):
        """
        Constructor.

        """
        self.align_input = alignfile
        self.target_fa   = targetfile
        self.query_fa    = queryfile
        self.opt         = options


class VariantsGenotype(object):
    """
    A class for variants' genotyping.

    Input required:
        (1) VCF files [list], bgzip format and has been index(.tbi) by tabix
        (2) Short reads' alignment files[list], sam/bam format and hash been
            sorted and index (.bai)
        (3) Reference fasta, has been index(.fai)
    
    """
    def __init__(self, vcffiles, bamfiles, referencefasta, options = None):
        """
        Constructor.

        Args:
            `vcffiles`: VCF files list.
            `bamfiles`: Bam files list. format: [(sample1,bam1), (), ...]
        """
        # Perhaps I should check the files before I open them?!
        # Checking whether vcffiles been tabix or not
        # Checking whether bamfiles been index or not
        # Checking whether fastafile been tabix or not
        self.ref_fasta   = pysam.FastaFile(referencefasta)
        self.vcfs        = vutil.VariantCandidateReader(vcffiles)
        # self.bam_readers is a hash dict, it's 'SampleID' point to
        # [index, bamReader] each bamfile represent to one sample
        self.bam_readers = {bf[0]:[i, pysam.AlignmentFile(bf[1])] 
                            for id, bf in enumerate(bamfiles)}
        # Sample's index => SampleID
        self.sample_index = {v[0]:k for k, v in self.bam_readers.items()}
        self.opt = options

    def genotyping(self):
        """
        The genotype process for haplotypes.
        """
        # loop all the vatiant windows
        for winvar in self._windowing():
            """
            Genotype the haplotypes in each windows and output the VCF
            """
            # Gernerate a list of haplotypes by all the combination of `winvar`
            haplotypes = gnt.generateAllHaplotypeByVariants(self.ref_fasta,
                                                            self.opt.max_read_len,
                                                            winvar)

            # Likelihood of each individual for all genotypes in the window.
            # The 'row' represent to genotypes  
            # The 'colum' represent to individuals
            # And it's a numpy array
            # CAUSION: These are variants' genotype likelihood, each value is
            # a likelihood for each individual in one diploid.
            # `genotype_likelihoods`: Recorde the likelihood of each 
            #                         genotype [It's not a log value now]
            # `genotype_hap_hash_id`: Recorde the genotypes by haplotypes' hash
            #                         id. [hap1_hash_id, hap2_hash_id]
            # `sample_map_nread`: It's used for recording the count of mapped
            #                     reads of this sample. And it's 2-d array 
            #                     [gnt][sample] 
            # See!! We don't have to reture the 'genotype' value, we just need
            # 'genotype_hap_hash_id' and 'haplotypes' , then we could get back
            # all the genotype easily! 
            genotype_likelihoods, genotype_hap_hash_id, sample_map_nread = (
                self._set_genotype_likelihood(haplotypes))

            # Now move to the next step. Use EM to calculate the haplotype
            # frequence in population scale.
            haplotypes_freq = self._calHaplotypeFreq(genotype_likelihoods,
                                                     genotype_hap_hash_id,
                                                     haplotypes)

            # convert to np.array.
            genotype_hap_hash_id = np.array(genotype_hap_hash_id)
            sample_map_nread     = np.array(sample_map_nread)

            # now we should start to pick the best genotype by finding the
            # max genotype of each individual in this window
            best_gnt_index       = genotype_likelihoods.argmax(axis = 0)
            individual_gnt_calls = genotype_hap_hash_id[best_gnt_index]
            # For the individual which have no reads aligning in this region
            individual_gnt_calls[sample_map_nread == 0] = None # No genotype
            
            # Now calcute the posterior probability for each 'winvar' variant
            # Type: dict(chrom = chrom, start = start, end = end, variant = var)
            # 'var_posterior_prob' is all the posterior porbability of variants
            # in `winvar['variants']`
            
            # I think var_posterior_prob shuld be a dict value. and the key is
            # the variant, the value is the posterior!
            var_posterior_prob = self.calVarPosteriorProb(haplotypes,
                                                          haplotypes_freq, 
                                                          genotype_hap_hash_id)

    def calVarPosteriorProb(self, haplotypes, hap_freq, genotype_hap_hash_id):
        """
        """
        posterior = {}
        done_var  = []
        for hap in haplotypes:
            for v in hap.variants:
                # I think var_posterior_prob shuld be a dict value. and the
                # key is the variant, the value is the posterior!
                post_prob    = self.calPosterior(v)
                posterior[v] = post_prob 

    def _calHaplotypeFreq(self, genotype_likelihoods, genotype_hap_hash_id,
                         haplotypes):
        """
        Calculate the haplotype frequence.
        """
        hap_num       = len(haplotypes)
        init_hap_freq = 1.0 / hap_num

        # Use haplotype's hash id to tracking the index of each haplotype
        hap_idx  = {hash(h):i for i, h in enumerate(haplotypes)}
        hap_freq = np.array([init_hap_freq for i in range(hap_num)])
        
        # Start to call EM algorithm.
        eps     = min(1e-3, 1.0 / (2 * hap_num))
        maxdiff = np.inf
        niter   = 0
        while maxdiff > eps and niter < DM.max_iter_num:
            # `hap_freq` will be updated automaticly when we call self.EM().
            maxdiff = self._EM(genotype_likelihoods, genotype_hap_hash_id,
                               hap_idx, hap_freq)
            niter  += 1

        return hap_freq
    
    def _EM(self, genotype_likelihoods, genotype_hap_hash_id, h_idx, hap_freq):
        """
        Perform one EM update. The update result will still store in `hap_freq`
        """
        # E Step: Estimate the expeceted values.
        emlikelihood = []
        for i, (h1, h2) in enumerate(genotype_hap_hash_id): # loop genotypes
            # Hardy-Weibery law to calculate "hp"
            # For 'homozygote'  : hp = hap1_freq * hap2_freq
            # For 'heterozygote': hp = 2 * hap1_freq * hap2_freq
            hp = hap_freq[h_idx[h1]] * hap_freq[h_idx[h2]] * (1 + (h1 != h2))

            tmp_lh = [glh * hp for glh in genotype_likelihoods[i]] 
            emlikelihood.append(tmp_lh)

        # Normalisation the genotype likelihood
        emlikelihood = self._normalisation(emlikelihood)

        # M Step: Re-estimate parameters.
        tmp_freq = np.zeros(len(hap_freq)) # Initial to 0.0
        for i, (h1, h2) in enumerate(genotype_hap_hash_id): # loop genotypes
            
            for lk in emlikelihood[i]:
                tmp_freq[h_idx[h1]] += lk
                tmp_freq[h_idx[h2]] += lk

        tmp_freq /= len(hap_freq) # From frequence to probability
        maxdiff   = np.abs(tmp_freq - hap_freq).max()
        hap_freq  = tmp_freq # OK, update the row haplotype frequence now!

        return maxdiff

    def _set_genotype_likelihood(self, haplotypes):
        """
        Setting the genotypes of the diploid combinated by 'haplotypes' 
        and return.

        Args:
            `haplotypes`: A list of Haplotype

        return a numpy array with individual likelihood for each genotype
        """
        genotypes = gnt.generateAllGenotypes(haplotypes)
 
        # Record buffer prevents to calculating again and saving time. 
        # But it must be a temporary value and clean the old data for each
        # genotype, or it'll cost a huge memory and it'll absolutely run  
        # out all the machine memory.  
        read_buffer_dict = {}
        # Record the genotype likelihood for genotyping 
        individual_loglikelihoods = []
        # Record the two haplotypes' hash id of each genotype
        genotype_hap_hash_id = []
        sample_map_nread     = []
        for gt in genotypes: 
            # `gt` is a list: [diploid, hap1_hash_id, hap2_hash_id]
    
            # Calculte the genotype likelihood for each individual
            individual_loglikelihoods, sample_nread = (
                self._calLikelihoodForIndividual(gt[0], read_buffer_dict))

            # The 'row' represent to each genotypes 
            # The 'colum' represent to each individuals 
            genotype_loglikelihoods.append(individual_loglikelihoods)
            genotype_hap_hash_id.append([gt[1], gt[2]])

            # read count each genotype/sample 
            if not sample_map_nread:
                # Assign one time is enough, they'll all be the same!
                sample_map_nread = sample_nread

        # Rescale genotype likelihood and covert to numpy array for the
        # next step. And causion: the array may contain value > 0 here,
        # before re-scale. They are all probability now after rescaling.
        # [NOW THEY ARE NOT log10 value any more!!] It's a 2D-array
        genotype_likelihoods = self._reScaleLikelihood(individual_loglikelihoods)

        return genotype_likelihoods, genotype_hap_hash_id, sample_map_nread
    
    def _normalisation(self, probability):
        """
        Nomalisation the probablity by using sum up the values of 
        each individual.
        """
        probability = np.array(probability, dtype = float)

        # Sum up the log likelihood for each individual of all genotypes
        sum_prob = probability.sum(axis = 0) # Sum up value by colums
        if sum_prob > 0:
            probability = probability / sum_prob # Normalisation 

        return probability

    def _reScaleLikelihood(self, likelihoods):
        """
        Re-scale the genotype log likelihoods by using maximum value
        """
        max_log = -1e7 # Threshold for maxmun likelihood value
        likelihoods = np.array(likelihoods, dtype = float)

        # We must assign the None to be 0.0. I choice 0.0 instead of other
        # value, cause I think we should just treat the probability to be 
        # 1.0 if the likelihood is None in individual of specify genotype.
        # [WE JUST DO IT HERE]
        lh_isnan = np.isnan(likelihoods)
        likelihoods[lh_isnan] = 0.0

        # Get the max log likelihood for each individual of all genotypes
        max_loglk = likelihoods.max(axis = 0) # Find the max value by colums
        # If all the likelihoods are small for some specify individuals, we
        # should still guarantee they're still small even after re-scaling.
        max_loglk[max_loglk < max_log] = max_log 

        # Re-scale by maxmum likelihood and the values are probability 
        return np.power(10, likelihoods - max_loglk)

    def _calLikelihoodForIndividual(self, genotype, read_buffer_dict):
        """
        Calculate the likelihood for each individual of this genotype by 
        re-aligning the reads of the specific individual.

        return an array loglikelihood for this genotype
        """

        # Initial likelihood to be None
        likelihoods = [None for i in self.sample_index]
        read_counts = [None for i in self.sample_index]
        for _, br in self.bam_readers.items():
            # Each bam file represent to one sample. Get the likelihood
            # by reads realignment

            # `lh` is log10 value and `rc` is the count of mapping reads
            lh, rc = genotype.calLikelihood(read_buffer_dict, br[1])
            # Storing by individual index. br[0] is the index for individual
            likelihoods[br[0]] = lh
            read_counts[br[0]] = rc

        # 'likelihood' is a 1D-array of log10 value
        return likelihoods, read_counts

    def _windowing(self):
        """
        Making windows along the reference genome and fill variants.
        Cutting the reference into windows by specify size, and fill it
        with variants

        return a list of variants 
        """
        varlist = []
        for chrom in self.ref_fasta.references:
            varlist += self._windowVarInSingleChrom(chrom)
        
        # returning all the variants in the VCF files
        return varlist

    def _windowVarInSingleChrom(self, chrom):    
        """
        Make windows for a specific chromosome of reference and fill 
        variants in each of them.

        Args:
            `chrom`: Chromososome id of reference fasta
        """
        if chrom not in self.ref_fasta.references:
            raise ValueError('"%s" not in the reference' % chrom)

        windows_varlist = []
        chrom_length = self.ref_fasta.get_reference_length(chrom)
        for start in range(1, chrom_length, self.opt.win_size): 

            end = min(start + self.opt.win_size - 1, chrom_length)
            var = self.vcfs.variants(chrom, start, end, self.opt.nosnp)

            # No variants in this window
            if not var: continue
           
            # Update the region's boundary by variant list 
            start = min([v.POS for v in var])
            end   = max([v.POS for v in var])

            # Yields a dictionary to store variants in this window. 
            wvar = dict(chrom = chrom, start = start, end = end, variant = var)
            if len(windows_varlist) == 0:
                windows_varlist.append(wvar)
                continue

            vn1 = len(var) # Variant number of this windows
            if vn1 > self.opt.max_var_num:
                # Too many variants in a window, spliting into samll pieces!
                for i in range(0, vn1, self.opt.max_var_num):
                    mi    = i + self.opt.max_var_num # Last index
                    start = min([v.POS for v in var[i:mi]])
                    end   = max([v.POS for v in var[i:mi]])
                    wvar  = dict(chrom = chrom, start = start, 
                                 end   = end, variant = var[i:mi])
                    windows_varlist.append(wvar)
            else: 

                vn2 = len(windows_varlist[-1]['variant']) # Variant's number
                # Try to merge window
                if start - windows_varlist[-1]['end'] <= self.opt.min_var_dist:
                    if vn1 + vn2 <= self.opt.max_var_num:
                        # Merging
                        windows_varlist[-1]['variant'] += var
                        windows_varlist[-1]['end']      = end
                    else:
                        windows_varlist.append(wvar)
                else:
                    windows_varlist.append(wvar)

        return windows_varlist


class VariantRecalibration(object):

    def __init__(self):
        """
        Constuctor.
        """
        self.vcfs = None

