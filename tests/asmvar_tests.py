from nose.tools import *
from pysam import FastaFile
import pysam

import asmvar.variantutil as vutil
from asmvar.variantutil import VariantCandidateReader as vcreader
from asmvar.haplotype import Haplotype as Hap
import asmvar.datum as dm
import asmvar.common as com
from asmvar.read import Read

def setup():
    print "SETUP!"

"""
def teardown():
    print "TEAR DOWN!"

def test_basic():
    print "I RAN!"
"""


def test_VariantCandidateReader():

    vcf_readers = vcreader(['tests/data/ex1.vcf'], 'options')
    vcf_readers = vcreader(['tests/data/tb.vcf.gz'], 'options')
    varlist   = vcf_readers.variants('20', 14369, 17330)
    test_case = [14370, 17330]
    series    = [i.POS for i in varlist]
    assert_equal(series, test_case)

    varlist = vcf_readers.variants('20')

def test_VariantCandidateReader_variants():

    vcf_readers = vcreader(['tests/data/ex1.vcf.gz'], 'options')
    varlist = []

    varlist = vcf_readers.variants('chr1')
    pos_case_chr1 = [288, 548, 1294]
    ref_case_chr1 = ['A', 'C', 'A']
    alt_case_chr1 = ['ACATAG', 'A', 'G']

    assert_equal(pos_case_chr1, [i.POS for i in varlist])
    assert_equal(ref_case_chr1, [i.REF for i in varlist])
    assert_equal(alt_case_chr1, [i.ALT[0].sequence for i in varlist])


    varlist = vcf_readers.variants('chr2')
    pos_case_chr2 = [156, 505, 784, 1344]
    ref_case_chr2 = ['A', 'A', 'C', 'A']
    alt_case_chr2 = ['AAG', 'G', 'CAATT', 'C']

    assert_equal(pos_case_chr2, [i.POS for i in varlist])
    assert_equal(ref_case_chr2, [i.REF for i in varlist])
    assert_equal(alt_case_chr2, [str(i.ALT[0]) for i in varlist])
    assert_equal(alt_case_chr2, [i.ALT[0].sequence for i in varlist])

    print '\n',str(i.ALT[0]) + 'A', '\t', i.ALT[0].sequence, '\n'
    
def test_vutil_get_sequence_context():

    fa = FastaFile('tests/data/ex1.fa')
    vcf_readers = vcreader(['tests/data/ex1.vcf.gz'], 'options')
    varlist = vcf_readers.variants('chr2')
    vutil.get_sequence_context(fa, varlist[0])

def test_vutil_homoRunForOneVariant():

    assert_equal(vutil._calHrunSize('tcggg'), 0)
    assert_equal(vutil._calHrunSize('ttcggg'), 2)
    assert_equal(vutil._calHrunSize('AATTGAGACTACAGAGCAAC'), 2)
    assert_equal(vutil._calHrunSize('ACTCACAGGTTTTATAAAAC'[::-1]), 0)

    fa = FastaFile('tests/data/ex1.fa')
    vcf_readers = vcreader(['tests/data/ex1.vcf.gz'], 'options')
    varlist = vcf_readers.variants('chr1')
    vutil.homoRunForOneVariant(fa, varlist[0])

    varlist = vcf_readers.variants('chr2')
    assert_equal(784, varlist[2].POS)
    assert_equal('ACTCACAGGTTTTATAAAAC', fa.fetch('chr2', varlist[2].POS - 20, varlist[2].POS))
    assert_equal('AATTGAGACTACAGAGCAAC', fa.fetch('chr2', varlist[2].POS, varlist[2].POS + 20))
    assert_equal('ACTCACAGGTTTTATAAAACAATTGAGACTACAGAGCAAC', fa.fetch('chr2', varlist[2].POS - 20, varlist[2].POS + 20))

    hr = vutil.homoRunForOneVariant(fa, varlist[2])
    assert_equal(2, hr)

    varlist[2].POS = fa.get_reference_length('chr2')
    hr = vutil.homoRunForOneVariant(fa, varlist[2])
    assert_equal(0, hr)

def test_datum():
    
    comdata = dm.CommonDatum()
    print 'indel_error_qual: ', comdata.homopol_penalty
    print 'hashmer: ', comdata.hashmer
    print 'hashsize: ', comdata.hashsize
    print 'max_align_size: ', comdata.max_align_size
    print 'indel_error_qual: ', comdata.homopol_penalty, '\n', dm.CommonDatum().homopol_penalty

def test_read():

    bam = pysam.AlignmentFile('tests/data/ex1.bam')
    read = Read()
    reads = []

    print '\n'
    for r in bam.fetch('chr1', 99, 100):
        reads.append(Read(r))
        print reads[-1].name, reads[-1].seqs, reads[-1].qual, reads[-1].seq_hash, len(reads[-1])
    


def test_common_SeqHashTable():

    ht  = com.SeqHashTable('') # Empty
    seq = 'ATCGCCGcccNatcgccgcccc'
    # Build the hash table for 'seq'
    ht  = com.SeqHashTable(seq, dm.CommonDatum().hashmer)

    print '\n', ht.hash_table, '\n'
    for id in ht.hash_pointer:
        print id, '=>', ht.hash_table[id]

    print '\n'
    idx = {}
    for id in ht.hash_pointer:
        idx[id] = idx.get(id, -1) + 1
        print idx[id], id, ht.hash_table[id][idx[id]]

def test_set_gap_open_cost(): 

    penalty = com.set_gap_open_penalty('ATCGCCGcccNatcgccgcccc', dm.CommonDatum().homopol_penalty)
    print '\npenalty:', penalty, [ord(i) for i in penalty]

    penalty = com.set_gap_open_penalty('atcg', dm.CommonDatum().homopol_penalty)
    print 'penalty:', penalty, [ord(i) for i in penalty]

    penalty = com.set_gap_open_penalty('AaAaaA', dm.CommonDatum().homopol_penalty)
    print 'penalty:', penalty, [ord(i) for i in penalty]

def test_Haplotype():

    fa  = FastaFile('tests/data/ex1.fa')
    hap = Hap(fa, 'chr1', 1, 20, 100)


def test_align_fastAlignmentRoutine():

    import ctypes
    ########
    class AlignTagPointer(ctypes.Structure):
        _fields_ = [("score", ctypes.c_int), ("pos", ctypes.c_int)]
    ########
    align = ctypes.CDLL('asmvar/align.so')

    print '\n\n** Now testing the align module **\n\n'

    seq1 = 'AAAGGGCAGGGGGGAGCACTAATGCGACCTCCACGCCCTTGTGTGTCCATGTACACACGCTGTCCTATGTACTTAT'
    #seq1 = 'AAAGGGCAGGGGGGAGCACTAATGCGACCTCCACGCCCTTGTGTGTGA'
    #seq2 = 'GGGAACAGGGGGGTGCACTAATGCGCTCCACGCC'
    seq2 = 'CAGGGGGGAGCACTAATGCGACCTCCACGCCCTTGT'
    qual = '<<86<<;<78<<<)<;4<67<;<;<74-7;,;8,;9'

    print seq1,'\n',seq2,'\n\n'

    aln1 = ''.join([str('\0') for i in range(2 * len(seq2) + 15)])
    aln2 = ''.join([str('\0') for i in range(2 * len(seq2) + 15)])
    local_gap_open = 'NKJHFA=854210/.-,,+**))(((\'\'\'&&&%%%$$$$#####"""""'
    
    read_start_in_hap = 0
    align.fastAlignmentRoutine.restype = ctypes.POINTER(AlignTagPointer)
    score = align.fastAlignmentRoutine(seq1, seq2, qual, len(seq2) + 15, len(seq2), 3, 2, local_gap_open, aln1, aln2)


    if score.contents.score == -1:
        score.contents.score = 1000
    print '\n\n*** After Align ***'
    print 'align score: ', score.contents.score, '\t', score.contents.pos
    print 'align1: ', aln1, len(aln1)
    print 'align2: ', aln2, len(aln2)
    print '\n'

    score.contents.score -= align.calculateFlankScore(len(seq1), 10, qual, local_gap_open, 3, 2, score.contents.pos + read_start_in_hap, aln1, aln2)

    #seq2 = 'CAGGGGGGAACTAATGCGACCTCCACGCCCTTAGTG'
    #score = align.fastAlignmentRoutine(seq1, seq2, qual, len(seq2) + 15, len(seq2), 3, 2, local_gap_open, aln1, aln2)
    print 'align score: ', score.contents.score, '\t', score.contents.pos
    print 'align1: ', aln1, len(aln1)
    print 'align2: ', aln2, len(aln2)
    print '\n'




