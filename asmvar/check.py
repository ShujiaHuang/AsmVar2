"""
This module use for checking the enviroment parameters
and checking the files: Fasta, bamfiles or vcf files, whether they 
could satisfy the requirement of AsmVar2.0 or not
"""
try:
    import pysam
except ImportError:
    raise Exception('"pysam" not available, try "pip install pysam"?')

try:
    import vcf
except ImportError:
    raise Exception('"vcf" not available, try "pip install PyVCF"?')

try:
    import numpy
except ImportError:
    raise Exception('"numpy" not available, try "pip install numpy"?')
