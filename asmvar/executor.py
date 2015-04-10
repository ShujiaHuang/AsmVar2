"""
The executor of AsmVar2 and it just could be called by ``Asmvar.py``
"""
import pysam
import vcf

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
        (1) VCF files [list], bgzip format and has been index(.tbi) by tabix
        (2) Short reads' alignment files[list], sam/bam format and hash been
            sorted and index (.bai)
        (3) Reference fasta, has been index(.fai)
    
    """
    def __init__(self, vcffiles, bamfiles, referencefasta, options = None):
        """
        Constructor.

        """
        # Perhaps I should check the files before I open them?!
        # Checking whether vcffiles been tabix or not
        # Checking whether bamfiles been index or not
        # Checking whether fastafile been tabix or not
        self.ref_fasta   = pysam.FastaFile(referencefasta)
        self.vcf_readers = vutil.VariantCandidateReader(vcffiles)
        self.bam_readers = [pysam.AlignmentFile(bf) for bf in bamfiles]
        self.opt         = options
        




