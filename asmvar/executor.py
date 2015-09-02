"""
This module will contain all the executor steps of variant calling,
genotying, variants' recalibration et.el

We have many important modules in AsmVar while this module is the lord
to rule them all, in a word, it's "The Ring". 
`Asmvar.py` is "Sauron", and this module could just be called by it.
"""
from genotype import genotyping as gnt

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
    A class for executing variants' genotyping.

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
        # executing genotyping for each variant in 'vcffile'
        exe_gnt = gnt.Genotype(vcffile, bamfiles, referencefasta, options)
	    # Calculat the variants' genotype and output
        exe_gnt.genotyping()


class VariantRecalibration(object):

    def __init__(self):
        """
        Constuctor.
        """
        self.vcfs = None

