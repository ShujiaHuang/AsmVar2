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
        # Index => SampleID
        self.sample_index = {v[0]:k for k, v in self.bam_readers.items()}
        self.opt = options

    def genotyping(self):
        """
        The genotype process for haplotypes.
        """
        # loop all the vatiant windows
        for winvar in self._windowing():
            """
            Genotype the haplotypes in each windows.
            """
            # Likelihood of each individual for all genotypes in the window.
            # The 'row' represent to each genotypes  
            # The 'colum' represent to each individuals
            # And it's a numpy array
            genotype_likelihoods = self.set_genotype_likelihood(winvar)

            # Now move to the next step. Use EM

    def set_genotype_likelihood(self, winvar):
        """
        Setting the genotypes for haplotyp of 'winvar' and return.

        Args:
            `winvar`: It's dict(chrom=chrom, start=start, end=end, variant=var)

        return a numpy array with individual likelihood for each genotype
        """
        genotypes = gnt.generateAllGenotypes(self.ref_fasta, 
                                             self.opt.max_read_len,
                                             winvar)
        # Record buffer prevents to calculating again and saving time. 
        # But it must be a temporary value and clean the old data for each
        # genotype, or it'll cost a huge memory and it'll absolutely run  
        # out all the machine memory  
        read_buffer_dict = {}
        # Record the genotype likelihood for genotyping 
        genotype_likelihoods = []
        for gt in genotypes: 
    
            # Calculte the genotype likelihood for each individual
            individual_loglikelihoods = self._calLikelihoodForIndividual(
                gt, read_buffer_dict) # A array of log10 value

            # The 'row' represent to each genotypes 
            # The 'colum' represent to each individuals 
            genotype_loglikelihoods.append(individual_loglikelihoods)

        # Rescale genotype likelihood and covert to numpy array for the
        # next step. And causion: the array may contain value > 0 here,
        # before re-scale. It's log10 value 
        genotype_likelihoods = self._reScaleLikelihood(genotype_likelihoods)

        return genotype_likelihoods

    def _reScaleLikelihood(self, likelihoods):
        """
        Re-scale the genotype log likelihoods by using maximum value
        """
        max_log = -1e6 # Threshold for maxmun likelihood value
        likelihoods = np.array(likelihoods, dtype = float)

        # Get the max log likelihood for each individual among all genotypes
        max_loglk = likelihoods.max(axis = 0) # Find max value by colums
        # If all the likelihoods are small for some specify individuals, we
        # should still guarantee they're still small even after re-scaling.
        max_loglk[max_loglk < max_log] = max_log 

        return likelihoods / max_loglk # Re-scale by maxmum likelihood

    def _calLikelihoodForIndividual(self, genotype, read_buffer_dict):
        """
        Calculate the likelihood for each individual of this genotype by 
        re-aligning the reads of the specific individual.

        return an array loglikelihood for this genotype
        """

        likelihoods = [None for i in self.sample_index] # initial array's size
        for _, br in self.bam_readers.items():
            # Each bam file represent to one sample. Get the likelihood
            # by reads realignment
            lh = genotype.calLikelihood(read_buffer_dict, br[1])
            # Storing by individual index. br[0] is the index for individual
            likelihoods[br[0]] = lh

        return likelihoods # A array of log10 value

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

