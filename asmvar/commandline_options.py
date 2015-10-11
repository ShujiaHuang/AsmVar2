"""
This is a user interface of AsmVar2. 
All the options and the usage for AsmVar2 will just be listed here
"""
import sys
import optparse

def genotype():
    """
    Loading all the parameters for genotyping.
    """

    usage = '\nUsage: python %prog genotype [options] > Output'
    optp = optparse.OptionParser(usage=usage)
    optp.add_option('-v', '--vcf', dest='vcffile', metavar='VCF', 
                    help = 'Variants. VCF format.', default='')
    optp.add_option('-b', '--bam', dest='bamfile', metavar='BAM', 
                    help = 'Bam Alignment file.', default='')
    optp.add_option('-r', '--ref', dest='ref_fasta_file', metavar='REF', 
                    help = 'Reference fa format.', default='')

    optp.add_option('-c', '--chr', dest='ref_chrom' , metavar='CHR', 
                    help = 'The chrom ID of Reference[ALL]. '
                           'eg: -c chr1 or -c chr1,chr2', default = '')
    optp.add_option('--nosnp', dest='nosnp', metavar='Boolen',
                    help = 'Remove all the SNPs in genotype', default=True)
    #optp.add_option('-w', '--win', dest='win_size', metavar='int', 
    #                help = 'Size of haplotype window', default=1000)
    optp.add_option('-d', '--variant_dist', dest='min_var_dist', metavar='int', 
                    help = 'The min variant distance windows', default=20)
    #optp.add_option('-n', '--num', dest='max_var_num', metavar='NUM',
    #                help = 'Max variants in window', default=1)
    optp.add_option('-s', '--read_len', dest='max_read_len', metavar='int',
                    help = 'Max length of reads', default=100)

    optp.add_option('-m', '--min_posterior', dest='min_posterior', metavar = 'int',
                    help = 'Minimum allowed value for somatic variant posterior', 
                    default = 5) 
    optp.add_option('-p', '--ped', dest='pedfile', metavar = 'STR', 
					help = 'pedigree information file. Not necessary.', 
                    default = '')
    optp.add_option('-f', '--fmt', dest='exfmt', metavar = 'STR',
					help = 'Add extra FORMAT fields from original VCF file. '
                           'And the parameter format could just be "AA:AC" '
                           'if you need it.', default = '')

    opt, _ = optp.parse_args()
    if not opt.vcffile: optp.error('Required [-v list of vcffile]\n')
    if not opt.bamfile: optp.error('Required [-b list of bamfile]\n')
    if not opt.ref_fasta_file: optp.error('Required [-r ref_fasta_file]\n')

    return opt
