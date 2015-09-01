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
import copy
import numpy as np

from multiprocessing import Pool 

import alignment as alg # The alignment module
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
    def __init__(self, vcffile, bamfiles, referencefasta, options = None):
        """
        Constructor.

        Args:
            `vcffile`: VCF file.
            `bamfiles`: Bam files list. format: [(sample1,bam1), (), ...]
        """
        # Perhaps I should check the files before I open them?!
        # Checking whether fastafile been tabix or not
        self.ref_fasta = pysam.FastaFile(referencefasta)
        self.vcfs      = vutil.VariantCandidateReader(vcffile)
        #self.bam_readers = {bf[0]:[i, pysam.AlignmentFile(bf[1])] 
        #                    for i, bf in enumerate(bamfiles)}

        # self.bam_readers is a hash dict, it's 'SampleID' point to
        # [index, bamfile] each bamfile represent to one sample
        self.sample2bamfile = {bf[0]: [i, bf[1]] for i, bf in enumerate(bamfiles)}
        # Sample's Index => SampleID
        self.index2sample = {v[0]: k for k, v in self.sample2bamfile.items()}
        self.samplenumber = len(self.index2sample)
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

        d = time.localtime()
        fileDate = '-'.join([str(d.tm_year), str(d.tm_mon), str(d.tm_mday)])
        vcf_header_info.record('##fileDate=' + fileDate)
        vcf_header_info.record('##reference=file://' + self.ref_fasta.filename)

        chrom   = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        samples = '\t'.join(s for k, s in self.index2sample.items())
        vcf_header_info.record(chrom + '\t' + samples)

        vcf_header_info.add('FORMAT', 'GT', 1, 'String', 'Genotype')
        vcf_header_info.add('FORMAT', 'GQ', 1, 'Integer', 'Genotype Quality')
        vcf_header_info.add('FORMAT', 'ND', '.', 'Integer', 'Number of reads '
                            'covering variant location in this sample')
        vcf_header_info.add('FORMAT', 'NV', '.', 'Integer', 'Number of reads '
                            'containing variant in this sample')
        vcf_header_info.add('FORMAT', 'PL', 'G', 'Integer',
            'Normalized, Phred-scaled likelihoods for genotypes as defined in '
            'the VCF specification for AA,AB and BB genotypes, where A = ref '
            'and B = variant. Only applicable for bi-allelic sites.')
        vcf_header_info.add('FORMAT', 'SB', 4, 'Integer', 
                            'Per-sample component statistics which comprise '
                            'the Fisher\'s Exact Test to detect strand bias. '
                            'The format: [ref-fwd, ref-reverse, non-ref-fwd, '
                            'non-ref-reverse]')
        vcf_header_info.add('FORMAT', 'AB', 1, 'Float', 'Allele balance for '
                            'each het genotype')

        vcf_header_info.add('INFO', 'HR', 1, 'Integer', 'Homozygous run')
        vcf_header_info.add('INFO', 'NR', 1, 'Float', 'N ratio around varant')
        return vcf_header_info

    def genotyping(self):
        """
        The genotype process for haplotypes.
        """
        # output vcf header information
        for k, h in sorted(self.vcf_header_info.header.items(), 
                           key = lambda d: d[0]): print h

        chr_id = self.ref_fasta.references # Defualt
        if self.opt.ref_chrom:
            chr_id = [c for c in self.opt.ref_chrom.split(',')]

        # loop all the vatiant windows
        for chrom in chr_id:

            if chrom not in self.ref_fasta.references:
                raise ValueError('"%s" not in the reference' % chrom)

            # Construct all the haplotype by VCF file
            # The format of hap_info's element is: 
            # [start, end, [haplotypes], wvar] 
            hap_info = self.construct_haplotype_by_variant(chrom)
            if not hap_info: continue
            # Calculate the likelihood information of haplotype in `hap_info`
            self.read_realign_to_haplotype(chrom, hap_info) 

            for h in hap_info:
                # h = [start, end, [haplotypes], wvar]
                # Genotyping each variant block and output
                self._exe_genotyping(h[2], h[3])

    def _exe_genotyping(self, haplotypes, winvar):
        """
        execute the genotyping process and output the result.

        Args:
            `haplotypes`: A list of haplotype
            `winvar`: All the variants in this haplotypes
        """
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
        genotypes = gnt.generateAllGenotypes(haplotypes)

        genotype_likelihoods = []
        genotype_hap_hash_id = []
        for gt in genotypes:
            genotype_likelihoods.append(gt.pop_loglikelihood)
            genotype_hap_hash_id.append([hash(gt.hap1), hash(gt.hap2)])

        genotype_likelihoods = self._reScaleLikelihoodInPop(genotype_likelihoods)
        # convert to np.array.
        genotype_hap_hash_id = np.array(genotype_hap_hash_id)
        sample_map_nread     = [len(lk) for lk in haplotypes[0].loglikelihood]
        sample_map_nread     = np.array(sample_map_nread)

        # Now move to the next step. Use EM to calculate the haplotype
        # frequence in population scale.
        hap_freq = self._calHaplotypeFreq(genotype_likelihoods,
                                          genotype_hap_hash_id,
                                          haplotypes)

        # now we should start to pick the best genotype by finding the
        # max genotype likelihood for each individual in this window
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

        """
        hap_idx  = {hash(h):i for i, h in enumerate(haplotypes)}
        for i, (h1, h2) in enumerate(genotype_hap_hash_id):# Debug
            h1 = hap_idx[h1]
            h2 = hap_idx[h2]
            print >> sys.stderr, "# Debug:", i, h1, h2, genotype_likelihoods[i], hap_freq[h1], hap_freq[h2]
        print >> sys.stderr, '\n'
        """
        #####################################################
        # All the prepare calculation are done! From now on #
        # we will calculate all the variant's posterior and #
        # output VCF                                        #
        #####################################################

        # Get getnotype and calculate the genotype posterior for 
        # each individuals on each positions
        hap_alt_var_hash = []
        hap_alt_var_cov  = {}
        for h in haplotypes:
            # The size of ALT in each haplotype's variants should just 
            # be one!
            tmp_hash_id = []
            for hv in h.variants:
                # The ALT could just have one variant in it!
                hs = hash((hv.CHROM, hv.POS, str(hv.ALT[0])))
                tmp_hash_id.append(hs)
                hap_alt_var_cov[hs] = hv.cov
            hap_alt_var_hash.append(tmp_hash_id)

        h_idx = {hash(h):i for i, h in enumerate(haplotypes)}
        for v in winvar:
            # Each variant means one position and remember there's REF
            # in v.ALT and REF is the first element.

            vcf_data_line        = vcfutils.Context()
            vcf_data_line.chrom  = v.CHROM
            vcf_data_line.pos    = v.POS
            vcf_data_line.Id     = v.ID
            vcf_data_line.ref    = v.REF
            vcf_data_line.alt    = [str(a) for a in v.ALT[1:]] # REF in ALT
            k = ':'.join([v.CHROM, str(v.POS), str(v.ALT[1])])
            vcf_data_line.qual   = var2posterior_phred[k]
            vcf_data_line.filter = '.'
            vcf_data_line.info   = {'HR': 'HR=' + str(v.hrun), 
                                    'NR': 'NR=' + str(v.nratio)}
            vcf_data_line.format = ['GT'] + sorted(['GQ', 'PL', 'ND', 'NV', 
                                                    'AB', 'SB'])

            # The first value in `alt_hash` is REF's hash id
            v_alt_hash = [hash((v.CHROM, v.POS, str(a))) for a in v.ALT]
            for i, s in self.index2sample.items(): # loop individuals
                # 'lh' are individual likelihoods for all the genotype.
                # 'non_ref_p': posterior of non reference call.
                # 'ref_p': posterior of reference call.
                lh, non_ref_p, ref_p = self.calGenotypeLikelihoodForIndividual(
                    v_alt_hash, 
                    {hs:d[i] for hs, d in hap_alt_var_cov.items()},
                    hap_alt_var_hash, 
                    h_idx, 
                    hap_freq,
                    genotype_likelihoods[:,i],
                    genotype_hap_hash_id)
                # Get the max likelihood and the corresponding allele index
                # We'll use alllele index to make the 'GT' in VCF
                max_lh, ale1, ale2, ND, NV, rf, rr, vf, vr = lh[lh[:,0].argmax()]
                ref_sb, var_sb = [int(rf), int(rr)], [int(vf), int(vr)]
                var_gnt_posterior  = max_lh / lh[:,0].sum()
                phred_ref_p        = self._phred(ref_p)
                phred_non_ref_p    = self._phred(non_ref_p)
                phred_var_gnt_posterior = self._phred(var_gnt_posterior)
                
                ale1, ale2 = int(ale1), int(ale2)
                if ale1 == -1: ale1 = '.' # '-1' represent fail
                if ale2 == -1: ale2 = '.' # '-1' represent fail
                GT = [ale1, '/', ale2]
                PL = [int(round(-10 * np.log10(x / max_lh))) 
                      for x in lh[:,0]]

                # Normalization the genotype with max likelihood
                if len(v.ALT) == 2: # There's REF in 'ALT', so we must use 2
                    # Bi-allele

                    # Don't make any call if the non-ref posterior is too low
                    if (phred_non_ref_p < self.opt.min_posterior and 
                        phred_ref_p < self.opt.min_posterior):
                        GT = ['.', '/', '.']
                    elif phred_non_ref_p < self.opt.min_posterior:
                        GT = [0, '/', 0] # Call ref  
                    else:
                        pass

                ab = '.'
                if ale1 != '.' and ale2 != '.' and ale1 != ale2 and ND > 0:
                    # Hete genotype
                    ab = round(float(NV) / ND , 3)

                sb = ','.join([str(d) for d in ref_sb] + 
                              [str(d) for d in var_sb])
                tmp = dict(GT = ''.join(str(gt) for gt in GT),
                           AB = str(ab),
                           GQ = str(phred_var_gnt_posterior),
                           PL = ','.join(str(pl) for pl in PL),
                           ND = str(int(ND)),
                           NV = str(int(NV)),
                           SB = sb)
                sample = ':'.join(tmp[k] for k in vcf_data_line.format)
                vcf_data_line.sample.append(sample)
            # Output VCF line
            vcf_data_line.print_context()

    def _reScaleLikelihoodInPop(self, loglikelihoods):
        """
        Rescale the genotype likelihood in population and covert to numpy array
        for the next step. 
        Causion: the array may contain value > 0 here before re-scale. They 
                 will all be probability now after rescaling.

        Args:
            `loglikelihoods`: It's 2-D array, row is genotyoe, colum is sample.
                              each element is likelihood for different individual.
        """
        loglikelihoods = np.array(loglikelihoods, dtype = float)

        # [WE JUST DO IT HERE]
        lh_isnan = np.isnan(loglikelihoods)
        loglikelihoods[lh_isnan] = -1e20

        # Normalization
        max_loglk = loglikelihoods.max(axis = 0)

        # If all the likelihoods are small for some specify individuals, we
        # should still guarantee they're still small even after re-scaling.
        max_loglk[max_loglk < -1e10] = 0.0

        # Re-scale by sum likelihood and the values are probability 
        lk = np.power(10, loglikelihoods - max_loglk)
        lk[lk < COMDM.min_float] = COMDM.min_float
        return lk

    def construct_haplotype_by_variant(self, chrom):
        """
        """
        hap_info  = []
        ref_fasta = self.ref_fasta.fetch(chrom)
        for wvar in self.load_vcf_by_chrom(ref_fasta, chrom):
            # Gernerate a list of haplotypes by all the combination of 
            # `wvar`.
            haplotypes = gnt.generateAllHaplotypeByVariants(
                ref_fasta, self.opt.max_read_len, wvar)

            # Pre-assign a 2D array to the `likelihood` for each hap
            # each row represent to sample, and each colum expected to
            # contain the read map likelihood of the specific sample
            # Use in `alg.alignReadToHaplotype`
            for h in haplotypes:
                h.loglikelihood = [[] for i in range(self.samplenumber)]

            # Construct all the haplotypes and the region in `haplotypes` 
            # list are all the same, they just dependant by `wvar`
            hap_info.append([haplotypes[0].hap_start, # hap region start 
                             haplotypes[0].hap_end,   # hap region end 
                             haplotypes,       # `haplotypes` is an list hap
                             wvar['variant']]) # All the variant in this hap

        return hap_info

    def read_realign_to_haplotype(self, chrom, hap_info):
        """
        Calculate the likelihood of each read map to each haplotype for
        each individual

        Args:
            `chrom`: the chromosome id
            `hap_info`: format is [start, end, [haplotypes]]
        """
        for i, s in self.index2sample.items(): 
            alg.alignReadToHaplotype(hap_info, 
                                     self.sample2bamfile[s][1], # The bamfile
                                     chrom, i)

        print >> sys.stderr, ('[INFO] %s Realignment and Haplotype likelihood '
                              'calculate done. Time: %s\n' % (chrom, time.asctime()))

    def load_vcf_by_chrom(self, chr_fa_seq, chrom):
        """
        Loading all the variant in the chromosome `chrom`
        """
        done_load_var = set()

        var_tmp_list = []
        for r in self.vcfs.vcf_reader.fetch(chrom):
            # Load variant each position
            if self.opt.nosnp and r.is_snp: continue

            h_id = hash((r.CHROM, r.POS, r.REF, len(r.ALT)))
            if (r.ALT[0] is None) or (h_id in done_load_var): continue
            done_load_var.add(h_id)

            for v in vutil.Variant(r).parse(): 
                v.hrun   = vutil.homoRunForOneVariant(chr_fa_seq, v) # Homo run
                v.nratio = vutil.nRatioForOneVariant(chr_fa_seq, v)
                var_tmp_list.append(v)
        var_tmp_list  = sorted(var_tmp_list) # Sorted by reference pos order
        variants_list = []

        for v in self._variantReConstruct(var_tmp_list):
            variants_list.append(v)

        print >> sys.stderr, '[INFO] Finish loading the variants of %s %s' % ( 
            chrom, time.asctime())
        # A list of dict, and the format for each element is: 
        # dict(chrom = chrom, start = start, end = end, variant = [var])
        return variants_list

    def _variantReConstruct(self, records):    
        """
        Args:
            records: A list of vcf.model._Record()
        """
        new_varlist = []

        for var in records:
            # Loop per variant
            start = var.POS
            end   = var.POS + len(var.REF) - 1

            # [2015-05-16] Put the reference sequence be the first element of ALT! 
            # This will be very convenient for us to use 0 => REF and other index
            # could be represent other variants. And this is very important in this
            # program.
            var.ALT = [vcf.model._Substitution(var.REF)] + var.ALT

            # Record coverage. Number of reads containing variant in different 
            # sample, the format is : [postive strand, reverse strand]
            var.cov = [[0, 0] for i in range(self.samplenumber)]

            # Yields a dictionary to store variants in this window. 
            wvar = dict(chrom = var.CHROM, start = start, end = end, variant = [var])
            new_varlist.append(wvar)

        return new_varlist # An array.

    def _phred(self, prob):
        """
        Calculate and return the phred score of 'prob'
        """
        phred = int(round(-10 * np.log10(max(1e-300, 1.0 - prob))))
        return min(200, phred)

    def calGenotypeLikelihoodForIndividual(self, 
                                           alt_hash_id, 
                                           hap_alt_var_cov_individual,
                                           hap_alt_var_hash, 
                                           h_idx, 
                                           hap_freq, 
                                           individual_genotype_likelihoods, 
                                           genotype_hap_hash_id):
        """
        Calculate genotype likelihoods for each individual.

        Args:
            `alt_hash_id`: The ALT hash id array, and the first element of ALT is REF!
            `hap_alt_var_cov_individual`: A dict. Variants' hash id to coverage!
            `hap_alt_var_hash`: A 2D array contain hash id of haplotype ALT
            `individual_genotype_likelihoods`: Specific individual's likelihood
                                               for all the genotypes.
            `genotype_hap_hash_id`: Hash id of two haplotypes in each genotype
        """
        non_ref_posterior, ref_posterior = 0.0, 0.0
        likelihoods    = []
        var_index      = range(len(alt_hash_id))
        individual_num = len(self.index2sample)

        ref_hash_id = alt_hash_id[0] # First element is REF's hash id
        # The index of 'var_index' could be use to represent REF or other 
        # variants, see in '_variantReConstruct'
        for i in var_index:
            var1_alt_hash = alt_hash_id[i]
            for j in var_index[i:]:
				# Marginal likelihood for this variant pair among all
				# the possible genotypes.
                var2_alt_hash = alt_hash_id[j]
                marginal_gnt_lh, nd, nv = 0.0, 0, 0
                ref_sb = [0, 0] # Record the reference strand bias
                var_sb = [0, 0] # Record the variant strand bias
				# Loop the different genotype by looping 'genotype_hap_hash_id'
                for k, (h1, h2) in enumerate(genotype_hap_hash_id):
                    
                    var1_in_hap1 = var1_alt_hash in hap_alt_var_hash[h_idx[h1]]
                    var1_in_hap2 = var1_alt_hash in hap_alt_var_hash[h_idx[h2]]
                    var2_in_hap1 = var2_alt_hash in hap_alt_var_hash[h_idx[h1]]
                    var2_in_hap2 = var2_alt_hash in hap_alt_var_hash[h_idx[h2]]

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
                    marginal_gnt_lh += current_lh

                    d1 = hap_alt_var_cov_individual[var1_alt_hash]
                    d2 = hap_alt_var_cov_individual[var2_alt_hash]

                    # Calculate ND and NV 
                    if var1_alt_hash != var2_alt_hash: 
                        nd  += (sum(d1) + sum(d2))
                    else: 
                        nd  += sum(d1)
                    
                    if var1_alt_hash != ref_hash_id:
                        nv += sum(d1)
                    if (var1_alt_hash != var2_alt_hash and 
                        var2_alt_hash != ref_hash_id): 
                        nv += sum(d2)

                    # Strand bias record
                    if var1_alt_hash == ref_hash_id:
                        ref_sb[0] += d1[0] # Forward 
                        ref_sb[1] += d1[1] # Reverse
                    else:
                        var_sb[0] += d1[0]
                        var_sb[1] += d1[1]

                    if var1_alt_hash != var2_alt_hash:
                        if var2_alt_hash == ref_hash_id:
                            ref_sb[0] += d2[0] # Forward
                            ref_sb[1] += d2[1] # Reverse
                        else:
                            var_sb[0] += d2[0]
                            var_sb[1] += d2[1]
                    ######## 
                    

                if alt_hash_id[j] != ref_hash_id or alt_hash_id[i] != ref_hash_id:
                    non_ref_posterior += marginal_gnt_lh
                else:
                    ref_posterior += marginal_gnt_lh

                ## Phased process
                phase_i, phase_j = -1, -1
                if alt_hash_id[j] == alt_hash_id[i]:
                    # Homo Ref or Homo Varirant. Do't need to phase
                    phase_i, phase_j = i, j
                elif alt_hash_id[j] == ref_hash_id or alt_hash_id[i] == ref_hash_id:
                    # Hete. Make sure call is phased correctly?
                    if var1_in_hap1:
                        phase_i, phase_j = i, j
                    elif var1_in_hap2:
                        phase_i, phase_j = j, i

                elif alt_hash_id[j] != alt_hash_id[i]:    
                    # Multi-allelic het. Make sure call is phased correctly?
                    if var1_in_hap1 and var2_in_hap2:
                        phase_i, phase_j = i, j 
                    elif var1_in_hap2 and var2_in_hap1:
                        phase_i, phase_j = j, i 
                # phase_i and phase_j could be used to represent the phased 
                # genotype.
                if marginal_gnt_lh < COMDM.min_float:
                    marginal_gnt_lh = COMDM.min_float 
                likelihoods.append([marginal_gnt_lh, phase_i, phase_j, 
                                    nd, nv, 
                                    ref_sb[0], ref_sb[1],  # Forward, Reverse 
                                    var_sb[0], var_sb[1]]) # Forward, Reverse

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
        sum_log10_prob = max(-300, np.log10(prob_mat.sum(axis = 0)).sum())

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

                k = ':'.join([v.CHROM, str(v.POS), str(v.ALT[0])])
                var2posterior_phred[k] = int(round(-10 * log_prob))

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


class VariantRecalibration(object):

    def __init__(self):
        """
        Constuctor.
        """
        self.vcfs = None

