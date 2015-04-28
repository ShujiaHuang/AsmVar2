"""
A class for record reads extracted from bamfile
"""
import pysam

class Read(object):

    def __init__(self, read = pysam.AlignedRead()):
        self.name = read.qname    # Read name
        self.seqs = read.query    # The read's sequence
        self.qual = read.qual     # The read's quality from bam alignment
        self.mapqual  = read.mapq # Mapping quality to reference genome by BWA
        self.seq_hash = None      # Hash the read's sequence
    
    def __len__(self):
        return len(self.seqs)

