"""
This module will contain all the executor steps of variant calling,
genotying, variants' recalibration et.el

We have many important modules in AsmVar while this module is the lord
to rule them all, in a word, it's "The Ring". 
`Asmvar.py` is "Sauron", and this module could just be called by it.
"""
import pysam
import vcf
import sys
import time
import numpy as np

import variantutil as vutil
import genotype as gnt
import vcfutils

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
                            for i, bf in enumerate(bamfiles)}
        # Sample's Index => SampleID
        self.index2sample = {v[0]:k for k, v in self.bam_readers.items()}
        self.opt          = options

        # Set vcf header
        self.vcf_header_info = self._set_vcf_header()

    def _set_vcf_header(self):
        """
        VCF Header infomation for output vcf file after genotyping
        """
        # Set the output vcf's header
        vcf_header_info = vcfutils.Header()
        vcf_header_info.record('##fileformat=VCFv4.1')

        d        = time.localtime()
        fileDate = '-'.join([str(d.tm_year), str(d.tm_mon), str(d.tm_mday)])
        vcf_header_info.record('##fileDate=' + fileDate)

        chrom   = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        samples = '\t'.join(s for k, s in self.index2sample.items())
        vcf_header_info.record(chrom + '\t' + samples)

        vcf_header_info.add('FORMAT', 'GT', 1, 'String', 'Genotype')
        vcf_header_info.add('FORMAT', 'GQ', 1, 'Integer', 'Genotype Quality')
        vcf_header_info.add('FORMAT', 'PL', 'G', 'Integer',
            'Normalized, Phred-scaled likelihoods for genotypes as '
            'defined in the VCF specification')

        return vcf_header_info

    def genotyping(self):
        """
        The genotype process for haplotypes.
        """
        # output vcf header information
        for k, h in sorted(self.vcf_header_info.header.items(), 
                           key = lambda d: d[0]): print h
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
            # `genotype_hap_hash_id`: Recorde haplotypes' hash id for each
            #                         genotypes. [hap1_hash_id, hap2_hash_id]. 
            #                         So 'genotype_hap_hash_id' is the same size
            #                         and order with genotype_likelihoods's colum
            # `sample_map_nread`: It's used for recording the count of mapped
            #                     reads of this sample. And it's 1-d array 
            #                     [sample], the sample order is the same with
            #                     the sample of `genotype_likelihoods`
            # See!! We don't have to reture the 'genotype' value, we just need
            # 'genotype_hap_hash_id' and 'haplotypes' , then we could get back
            # all the genotype easily! 
            genotype_likelihoods, genotype_hap_hash_id, sample_map_nread = (
                self._set_genotype_likelihood(haplotypes))

            # Now move to the next step. Use EM to calculate the haplotype
            # frequence in population scale.
            hap_freq = self._calHaplotypeFreq(genotype_likelihoods,
                                              genotype_hap_hash_id,
                                              haplotypes)

            # convert to np.array.
            genotype_hap_hash_id = np.array(genotype_hap_hash_id)
            sample_map_nread     = np.array(sample_map_nread)

            # now we should start to pick the best genotype by finding the
            # max genotype for each individual in this window
            best_gnt_index       = genotype_likelihoods.argmax(axis = 0)
            individual_gnt_calls = genotype_hap_hash_id[best_gnt_index]
            # For the individual which have no reads aligning in this region
            individual_gnt_calls[sample_map_nread == 0] = [-1, -1] # No genotype
            
            # Now calcute the posterior probability for each 'winvar' variant
            # Type: dict(chrom = chrom, start = start, end = end, variant = var)
            # 'var_posterior_prob' is all the posterior porbability of variants
            # in `winvar['variants']`
            
            # `var2posterior_phred` is a dict value. and the key is
            # the variant, the value is the posterior!
            var2posterior_phred = self.calVarPosteriorProb(haplotypes,
                                                           hap_freq, 
                                                           genotype_likelihoods,
                                                           genotype_hap_hash_id,
                                                           sample_map_nread)
            
            #####################################################
            # All the prepare calculation are done! From now on #
            # we will calculate all the variant's posterior and #
            # output VCF                                        #
            #####################################################

            # Get getnotype and calculate the genotype posterior for 
            # each individuals on each positions
            hap_alt_var_hash = []
            for h in haplotypes:
                # The size of ALT in each haplotype's variants should just 
                # be one!
                hap_alt_var_hash.append([hash((hv.CHROM,hv.POS,str(hv.ALT[0]))) 
                                         for hv in h.variants])

            h_idx = {hash(h):i for i, h in enumerate(haplotypes)}
            for v in winvar['variant']:
                # Each variant means one position and remember there's REF
                # in v.ALT and REF is the first element.

                vcf_data_line = vcfutils.Context()
                vcf_data_line.chrom = v.CHROM
                vcf_data_line.pos   = v.POS
                vcf_data_line.Id    = v.ID
                vcf_data_line.ref   = v.REF
                vcf_data_line.alt   = [str(a) for a in v.ALT[1:]] # REF in ALT
                vcf_data_line.qual  = v.QUAL
                vcf_data_line.filter = v.FILTER
                vcf_data_line.info   = v.INFO
                vcf_data_line.format = ['GT'] + sorted(['GQ', 'PL'])

                for i, s in self.index2sample.items(): # loop individuals
                    # 'lh' are individual likelihoods for all the genotype.
                    # 'non_ref_p': posterior of non reference call.
                    # 'ref_p': posterior of reference call.
                    lh, non_ref_p, ref_p = self.calGenotypeLikelihoodForIndividual(
                        v, hap_alt_var_hash, h_idx, hap_freq,
                        genotype_likelihoods[:,i],
                        genotype_hap_hash_id,
                        sample_map_nread[i])

                    # Get the max likelihood and the corresponding allele index
                    # We'll use alllele index to make the 'GT' in VCF
                    max_lh, ale1, ale2 = lh[lh[:,0].argmax()]
                    var_gnt_posterior  = max_lh / lh[:,0].sum()
                    phred_ref_p        = self._phred(ref_p)
                    phred_non_ref_p    = self._phred(non_ref_p)
                    phred_var_gnt_posterior = self._phred(var_gnt_posterior)
                    
                    ale1, ale2 = int(ale1), int(ale2)
                    if ale1 == -1: ale1 = '.' # '-1' represent fail phased
                    if ale2 == -1: ale2 = '.' # '-1' represent fail phased
                    GT = [str(ale1), '/', str(ale2)]
                    PL = [int(round(-10 * np.log10(max(x / max_lh, 1e-300)))) 
                          for x in lh[:,0]]
                    # Normalization the genotype with max likelihood
                    if len(v.ALT) == 2: # There's REF in 'ALT', so we must use 2
                        # Bi-allele
                        normarl_GLs = [round(np.log10(max(x / max_lh, 1e-300)), 2) 
                                       for x in lh[:,0]]

                        if (phred_non_ref_p < self.opt.min_posterior and 
                            phred_ref_p < self.opt.min_posterior):
                            # Don't make any call if the non-ref posterior 
                            # is too low
                            GT = ['.', '/', '.']
                        elif phred_non_ref_p < self.opt.min_posterior:
                            GT = ['0', '/', '0'] # Call ref  
                        else:
                            pass
                    else:
                        # multiple alleles
                        normarl_GLs = [-1, -1, -1]
                    tmp = dict(GT = ''.join(GT),
                               GQ = str(phred_var_gnt_posterior),
                               PL = ','.join(str(l) for l in PL))
                    sample = ':'.join(str(tmp[k]) for k in vcf_data_line.format)
                    vcf_data_line.sample.append(sample)
                # Output VCF line
                vcf_data_line.print_context()

    def _phred(self, prob):
        """
        Calculate and return the phred score of 'prob'
        """
        phred = int(round(-10 * np.log10(max(1e-300, 1.0 - prob))))
        return min(100, phred)

    def calGenotypeLikelihoodForIndividual(self, 
                                           variant, 
                                           hap_alt_var_hash, 
                                           h_idx, 
                                           hap_freq, 
                                           individual_genotype_likelihoods, 
                                           genotype_hap_hash_id,
                                           map_read_num):
        """
        Calculate genotype likelihoods for each individual.

        Args:
            `variant`: vcf.model._Record, and the first element of ALT is REF!
            `hap_alt_var_hash`: A 2D array contain hash id of haplotype ALT
            `individual_genotype_likelihoods`: Specific individual's likelihood
                                               for all the genotypes.
            `genotype_hap_hash_id`: Hash id of two haplotypes in each genotype
            `map_read_num`: Mapping reads number of this individual
        """
        alt_hash = [hash((variant.CHROM, variant.POS, str(a))) 
                    for a in variant.ALT]

        non_ref_posterior, ref_posterior = 0.0, 0.0
        likelihoods    = []
        var_index      = range(len(variant.ALT))
        individual_num = len(self.index2sample)

        # The index of 'var_index' could be use to represent REF or other 
        # variants, see in '_windowVarInSingleChrom'
        for i in var_index:
            for j in var_index[i:]:
				# Marginal likelihood for this variant pair among all
				# the possible genotypes.
                marginal_gnt_lh = 0.0
				# Loop the genotype by looping 'genotype_hap_hash_id'
                for k, (h1, h2) in enumerate(genotype_hap_hash_id):
                    
                    var1_in_hap1 = alt_hash[i] in hap_alt_var_hash[h_idx[h1]]
                    var1_in_hap2 = alt_hash[i] in hap_alt_var_hash[h_idx[h2]]
                    var2_in_hap1 = alt_hash[j] in hap_alt_var_hash[h_idx[h1]]
                    var2_in_hap2 = alt_hash[j] in hap_alt_var_hash[h_idx[h2]]

                    # Just accumulate likelihood for the satisfy genotype.
                    if ((not (var1_in_hap1 and var2_in_hap2)) and 
                        (not (var1_in_hap2 and var2_in_hap1))): continue

                    # Only use EM frequencies for large-ish populations.
                    if individual_num > 25:
                        hp = ((1 + (h1 != h2)) * hap_freq[h_idx[h1]] * 
                               hap_freq[h_idx[h2]])
                    else:
                        hp = (1 + (h1 != h2))

                    current_lh = hp * individual_genotype_likelihoods[k]
                    #print 'current_lh:', k, i, j, individual_genotype_likelihoods[k], current_lh
                    marginal_gnt_lh += current_lh

                if variant.ALT[j] != variant.REF or variant.ALT[i] != variant.REF:
                    non_ref_posterior += marginal_gnt_lh
                else:
                    ref_posterior += marginal_gnt_lh

                ## Phased process
                phase_i, phase_j = -1, -1
                if variant.ALT[j] == variant.ALT[i]:
                    # Homo Ref or Homo Varirant. Do't need to phase
                    phase_i, phase_j = i, j
                elif variant.ALT[j] == variant.REF or variant.ALT[i] == variant.REF:
                    # Het. Make sure call is phased correctly?
                    if var1_in_hap1:
                        phase_i, phase_j = i, j
                    elif var1_in_hap2:
                        phase_i, phase_j = j, i

                elif variant.ALT[j] != variant.ALT[i]:    
                    # Multi-allelic het. Make sure call is phased correctly?
                    if var1_in_hap1 and var2_in_hap2:
                        phase_i, phase_j = i, j 
                    elif var1_in_hap2 and var2_in_hap1:
                        phase_i, phase_j = j, i 
                # phase_i and phase_j could be used to represent the phased 
                # genotype.
                likelihoods.append([marginal_gnt_lh, phase_i, phase_j])
                #print '** marginal_gnt_lh:', phase_i, phase_j, marginal_gnt_lh, '\n'

        likelihoods        = np.array(likelihoods)
        non_ref_posterior /= likelihoods[:,0].sum()
        ref_posterior     /= likelihoods[:,0].sum()
        # likelihoods the genotypes of all the ALTs' combination on this
        # position of this specific individual
        return likelihoods, non_ref_posterior, ref_posterior

    def calVarPosteriorProb(self, haplotypes, hap_freq, genotype_likelihoods,
                            genotype_hap_hash_id, sample_map_nread):
        """
        Posterior for each variants
        """
        # Use haplotype's hash id to tracking the index of each haplotype  
        h_idx = {hash(h):i for i, h in enumerate(haplotypes)}
        prob_mat = self._cal_prob_mat(h_idx, hap_freq, genotype_likelihoods,
                                      genotype_hap_hash_id)

        # Assign the probability to be 0 if no reads covert the individual
        # And it's still a 2-D array.
        prob_mat *= (sample_map_nread > 0)
        prob_mat[prob_mat <= 0] = COMDM.min_float

        # We sum up all the probability of genotype for each individual, 
        # and then times all the sum-probability together.
        # It's a log10 value now.
        sum_log10_prob = np.log10(prob_mat.sum(axis = 0)).sum()

        var2posterior_phred = {}
        for hap in haplotypes:

            for v in hap.variants:

                if v in var2posterior_phred: continue
                # Calculate the new haplotype frequence if they do not have 
                # 'v' in it.
                hap_freq_without_v = self.cal_hap_freq_without_var(v,
                                                                   haplotypes,
                                                                   hap_freq)
                prob_mat_without_v = self._cal_prob_mat(h_idx, 
                                                        hap_freq_without_v,
                                                        genotype_likelihoods,
                                                        genotype_hap_hash_id)

                prob_mat_without_v *= (sample_map_nread > 0)
                prob_mat_without_v[prob_mat_without_v <= 0] = COMDM.min_float
                # It's a log10 value without 'var'
                sum_log10_prob_without_var = np.log10(
                    prob_mat_without_v.sum(axis = 0).sum())

                delta = max(-300, sum_log10_prob_without_var - sum_log10_prob)
                ratio = max(COMDM.min_float, 10 ** delta)

                prior = vutil.calPrior(self.ref_fasta, v)
                # The key is the variant, the value is the posterior_phred!
                # Remember: 'hap.variants' will been change if we change 'v'
                # The posterior probability is: 
                #   prior / (prior + ratio * (1.0 - prior))
                log_prob = (np.log10(ratio) + np.log10(1.0 - prior) - 
                            np.log10(prior + ratio * (1.0 - prior)))

                var2posterior_phred[v] = round(-10 * log_prob)

        # Return a dict record the variants' posterior phred score
        return var2posterior_phred # For each variants
    
    def cal_hap_freq_without_var(self, variant, haplotypes, hap_freq):
        """
        """
        hap_freq_without_var = np.zeros(len(hap_freq))
        sum_freq = 0.0 

        for i, h in enumerate(haplotypes):
    
            if variant not in h.variants:
                # Record all the hap_freq which not have 'variant' in it.
                hap_freq_without_var[i]  = hap_freq[i]
                sum_freq                += hap_freq[i]

        if sum_freq > 0:
            hap_freq_without_var /= sum_freq
        
        return hap_freq_without_var
        
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
        while maxdiff > eps and niter < COMDM.max_iter_num:
            # `hap_freq` will be updated automaticly after we call self.EM().
            niter += 1
            maxdiff, hap_freq = self._EM(genotype_likelihoods, 
                                         genotype_hap_hash_id,
                                         hap_idx, hap_freq)
        return hap_freq
    
    def _EM(self, genotype_likelihoods, genotype_hap_hash_id, h_idx, hap_freq):
        """
        Perform one EM update.
        """
        # E Step: Estimate the expeceted values.
        emlikelihood = self._cal_prob_mat(h_idx, hap_freq, 
                                          genotype_likelihoods,
                                          genotype_hap_hash_id)
        # Normalisation the genotype likelihood
        emlikelihood = self._normalisation(emlikelihood)

        # M Step: Re-estimate parameters.
        new_freq = np.zeros(len(hap_freq)) # Initial to 0.0
        for i, (h1, h2) in enumerate(genotype_hap_hash_id): # loop genotypes
            
            for lk in emlikelihood[i]:
                new_freq[h_idx[h1]] += lk
                new_freq[h_idx[h2]] += lk

        #new_freq /= len(hap_freq) # From frequence number to probability
        sample_num = max(1, len(genotype_likelihoods[0]))
        new_freq  /= 2 * sample_num
        maxdiff    = np.abs(new_freq - hap_freq).max()
        return maxdiff, new_freq

    def _cal_prob_mat(self, h_idx, hap_freq, genotype_likelihoods,
                      genotype_hap_hash_id):
        """
        """
        prob = [] 
        # Use haplotype's hash id to tracking the index of each haplotype  
        for i, (h1, h2) in enumerate(genotype_hap_hash_id): # loop genotypes
            # Hardy-Weibery law to calculate "hp"
            # For 'homozygote'  : hp = hap1_freq * hap2_freq
            # For 'heterozygote': hp = 2 * hap1_freq * hap2_freq
            hp = ((1 + (h1 != h2)) * hap_freq[h_idx[h1]] * hap_freq[h_idx[h2]])
            tmp_lh = [lh * hp for lh in genotype_likelihoods[i]]
            prob.append(tmp_lh)

        return np.array(prob) # It's a 2d array

    def _normalisation(self, probability):
        """
        Nomalisation the probablity by using the sum of the values of 
        each individual.
        """
        probability = np.array(probability, dtype = float)
        # Sum up the log likelihood for each individual of all genotypes
        sum_prob = probability.sum(axis = 0) # Sum up value by colums
        sum_prob[sum_prob == 0] = 1.0
        probability = probability / sum_prob # Normalisation 

        return probability

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
        genotype_loglikelihoods = []
        # Record the two haplotypes' hash id of each genotype
        genotype_hap_hash_id = []
        sample_map_nread     = []
        for gt in genotypes: 
            # `gt` is a list: [diploid, hap1_hash_id, hap2_hash_id]
            # Calculte the genotype likelihood for each individual
            # This is the most costly performance in this program!!
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
        # next step. And causion: the array may contain value > 0 here
        # before re-scale. They will all be probability now after rescaling.
        # [NOW THEY ARE NOT log10 value any more!!] 2D-array
        genotype_likelihoods = self._reScaleLikelihood(genotype_loglikelihoods)
        #for g in genotype_likelihoods: print 'genotype_likelihoods:', g

        return genotype_likelihoods, genotype_hap_hash_id, sample_map_nread
    
    def _reScaleLikelihood(self, likelihoods):
        """
        Re-scale the genotype log likelihoods by using maximum value
        """
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
        # [need this??] max_loglk[max_loglk < COMDM.min_log] = COMDM.min_log 

        # Re-scale by maxmum likelihood and the values are probability 
        lk = np.power(10, likelihoods - max_loglk)
        lk[lk < COMDM.min_float] = COMDM.min_float
        return lk

    def _calLikelihoodForIndividual(self, genotype, read_buffer_dict):
        """
        Calculate the likelihood for each individual of this genotype
        by re-aligning the reads of the specific individual.

        return an array loglikelihood for this genotype
        """
        # Initial likelihood to be None
        likelihoods = [None for i in self.index2sample]
        read_counts = [None for i in self.index2sample]
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
        chr_id = self.ref_fasta.references # Defualt
        if self.opt.ref_chrom:
            chr_id = [c for c in self.opt.ref_chrom.split(',')]

        varlist = []
        for chrom in chr_id:
            for v in self._windowVarInSingleChrom(chrom):
                varlist.append(v)
        
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

        done_load_var   = set() # A set record the loaded variants
        windows_varlist = []
        chrom_length = self.ref_fasta.get_reference_length(chrom)
        for start in range(1, chrom_length, self.opt.win_size): 

            var = []
            end = min(start + self.opt.win_size - 1, chrom_length)
            var = self.vcfs.variants(chrom, start, end, 
                                     self.opt.nosnp, done_load_var)
            # No variants in this window
            if not var: continue
            # Update the region's boundary by using the new variants' data 
            start = min([v.POS for v in var])
            end   = max([v.POS + len(v.REF) - 1 for v in var])

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
                    end   = max([v.POS + len(v.REF) - 1 for v in var[i:mi]])
                    wvar  = dict(chrom = chrom, start = start, 
                                 end   = end, variant = var[i:mi])
                    windows_varlist.append(wvar)
            else: 

                vn2 = len(windows_varlist[-1]['variant']) # Variant's number
                # Try to merge window
                if start - windows_varlist[-1]['end'] <= self.opt.min_var_dist:

                    if start < windows_varlist[-1]['end']:
                        # God! They're Overlap, do NOT merge !!
                        windows_varlist.append(wvar) 
                    elif vn1 + vn2 <= self.opt.max_var_num:
                        # Merging
                        windows_varlist[-1]['variant'] += var
                        windows_varlist[-1]['end']      = end
                    else:
                        windows_varlist.append(wvar)
                else:
                    windows_varlist.append(wvar)

            del_max_pos = max([v.POS + len(v.REF) - 1 for v in windows_varlist[-1]['variant']])

        # [2015-05-16] Put the reference sequence be the first element of ALT! 
        # This will be very convenient for us to use 0 => REF and other index
        # could be represent other variants. And this is very important in case
        # of making error.
        for i, wv in enumerate(windows_varlist):
            for v in wv['variant']:
                # Because the feature of Python, we can just change 'v', but
                # it will still record the change to windows_varlist
                v.ALT = [vcf.model._Substitution(v.REF)] + v.ALT

        return windows_varlist # An array.


class VariantRecalibration(object):

    def __init__(self):
        """
        Constuctor.
        """
        self.vcfs = None

