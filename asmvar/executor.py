"""
The executor of AsmVar2 and it just could be called by ``Asmvar.py``
"""
import variantutil as vutil
import haplotype

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
    


class VariantGenotype(object):
    """
    A class for variants' genotyping.

    Input required:
        (1) VCF files, bgzip format and has been index(.tbi) by tabix
        (2) Short reads' alignment files, sam/bam format and hash been
            sorted and index (.bai)
        (3) Reference fasta, has been index(.fai)
    
    """
    def __init__(self, vcffile, bamfile, referencefasta, options):
        """
        Constructor.

        """
        self.vcf_infile = vcffile
        self.bam_input  = bamfile
        self.ref_fasta  = referencefasta
        self.opt = options
        

