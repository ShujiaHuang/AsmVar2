"""
This module contain the alignment functions for reads to haplotype,
genotyping or realignment process
"""
import os
import ctypes
try:
    import pysam
except ImportError:
    pysam = None
    raise Exception('pysam not available, try "pip install pysam"?')

import common as com # Same common functions

# The diretory of this program
dir    = os.path.dirname(os.path.abspath(__file__))
align  = ctypes.CDLL(dir + '/align.so')  # The alignment module written by C

def alignReadToHaplotype(haplotype, bam_stream):
    """ Learn from Platypus.

    This is the most basic and important part of AsmVar2. This function decides
    where to anchor the read to the specified haplotype, and calls the banded-
    alignment function -- fastAlignmentRoutine function in 'align.c'. 

    Args:

    Return:

    Requires pysam
    """

    if haplotype.seq_hash is None:
        haplotype.seq_hash = com.SeqHashTable(haplotype.sequence)


    
def singleRead2Haplotype(haplotype, read):
    """
    Mapped a single read to a specified haplotype sequence. 

    This is used to find the best anchor position in the haplotype sequence
    for the specified read.

    Args:

    Return: The best alignement score
    """
    
    











