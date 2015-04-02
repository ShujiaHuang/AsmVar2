from nose.tools import *
from pysam import FastaFile

import asmvar.variantutil as vutil
from asmvar.variantutil import VariantCandidateReader as vcreader
from asmvar.haplotype import Haplotype as Hap

def setup():
    print "SETUP!"

def teardown():
    print "TEAR DOWN!"

def test_basic():
    print "I RAN!"

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
    assert_equal(alt_case_chr2, [i.ALT[0].sequence for i in varlist])
    
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

def test_Haplotype():

    fa  = FastaFile('tests/data/ex1.fa')
    hap = Hap(fa, 'chr1', 1, 20)

