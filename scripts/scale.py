"""
This is a program use for scale VCF file into small ones by 
custom parameter. And compress all the sub-file by bgzip and
build the tabix index.

Copyright (c) Shujia Huang
Date: 2015-08-11
"""
import sys
import os
import optparse
import string
from subprocess import Popen, PIPE

import vcf
import pysam
#from pysam import TabixFile # I use the tabix module in pysam

def get_opt():
    """ Loading parameter for scaling vcf. """
    usage = '\nUsage: python '
    optp = optparse.OptionParser(usage=usage)
    optp.add_option('-v', '--vcf', dest = 'vcffile', metavar = 'VCF', 
                    help = 'Variants. VCF format.', default = '')
    optp.add_option('-c', '--chr', dest='ref_chrom' , metavar='CHR',
                    help = 'The chrom ID of Reference[ALL]. e.g. '
                    '-c chr1 or -c chr1,chr2', default = '')
    optp.add_option('-n', '--num', dest = 'number', metavar = 'INT', 
                    help = 'The number you intend to split input by -v.',
                    default = 10)
    optp.add_option('-d', '--outdir', dest = 'outdir', metavar = 'DIR',
                    help = 'The subfile output directory.', default = '.')
    optp.add_option('-p', '--prefix', dest = 'outprefix', metavar = 'OUT',
                    help = 'Out file prefix', default = 'test')
    optp.add_option('--rec', dest = 'recursive', metavar = 'Boolen',
                    help = 'Recursvie to split the input vcffile(-v) by '
                    'its chrom ID(-c) or not.', default = '')

    opt, _ = optp.parse_args()
    if not opt.vcffile: optp.error('Required [-v vcffile]\n')
    print >> sys.stderr, 'Parameters: python', ' '.join(sys.argv) 

    opt.number = abs(string.atoi(opt.number))
    opt.recursive = True if opt.recursive else False
    if not os.path.exists(opt.outdir):
        os.makedirs(opt.outdir)

    return opt

def main(opt):

    if opt.number > 0 or opt.recursive:
        sub_vcf_files = splitVCF(opt.vcffile, opt.ref_chrom, opt.number,
                                 opt.outdir + '/tmp_in_dir', opt.recursive)
    else:
        # Don't need to split the VCF file
        sub_vcf_files = [opt.vcffile]

    print '\n'.join(sub_vcf_files)

def createJobScript(input_files):
    pass

def splitVCF(vcffile, ref_chrom, split_num, sub_outdir, is_rec_split = True):

    vcf_line_count, chrom_ids = _get_vcf_line_count(vcffile, ref_chrom)

    # get vcffile's file name by os.path.split
    _, fname = os.path.split(vcffile)
    if not os.path.exists(sub_outdir):
        os.makedirs(sub_outdir)

    vcf_reader = vcf.Reader(filename = vcffile)
    sub_vcf_files = []
    if is_rec_split: 
        """
        Split the whole vcf file by different chrom in `chrom_ids`.
        """
        for chrom in chrom_ids:

            sub_chr_dir = sub_outdir + '/' + chrom
            if not os.path.exists(sub_chr_dir):
                os.makedirs(sub_chr_dir)
                
            tot_num, lniof = _set_step_num(vcf_line_count[chrom], split_num)
            outprefix = sub_chr_dir + '/tmp.in.' + fname + '.' + chrom
            for f in outputSubVCF(vcf_reader.fetch(chrom), lniof, 
                                  tot_num, outprefix):
                sub_vcf_files.append(f)
    else:
        outprefix = sub_outdir + '/tmp.in.' + fname
        tot_num, lniof = _set_step_num(vcf_line_count['all'], split_num)
        sub_vcf_files  = outputSubVCF(vcf_reader, lniof, tot_num, outprefix)

    # compress all the sub file for tabix
    for i, f in enumerate(sub_vcf_files): 

        # Build the tabix index. The method will automatically compressed the 
        # file if which name does not end with '.gz' and the original file will
        # be removed and only the compressed file will be retained
        f = pysam.tabix_index(sub_vcf_files[i], force = True, preset = 'vcf')
        sub_vcf_files[i] = f # The compressed file after tabix

    return sub_vcf_files

def outputSubVCF(vcf_reader, line_num_in_one_file, t_f_n, 
                 outfile_prefix):
    """
    Split the whole vcf file into several sub-files. Compress by bgzip
    and creat the tabix index.

    return all the path of subfiles as a list.

    Args:
        `t_f_n`: Total file numbers. ('split_num' same as opt.number)
    """
    line_num = 0
    sfn      = 0
    sub_files  = []
    vcf_writer = None
    for r in vcf_reader:

        if line_num % line_num_in_one_file == 0:

            if vcf_writer:
                vcf_writer.close()

            sfn += 1
            sub_file = outfile_prefix + '.' + str(sfn) + '_' + str(t_f_n) + '.vcf'
            vcf_writer = vcf.Writer(open(sub_file, 'w'), vcf_reader)
            sub_files.append(sub_file)

        # Output the VCF record
        vcf_writer.write_record(r)

        line_num += 1

    if vcf_writer:
        vcf_writer.close()

    return sub_files

def _get_vcf_line_count(vcffile, chrom_id):
    """
    Get line number of this vcf file. And record them in a dict.
    """

    f = pysam.TabixFile(vcffile)
    chrom_id_set = set(f.contigs)
    f.close()

    if chrom_id:
        chrom_ids = set(chrom_id.split(','))

    line_count = {}
    vcf_reader = vcf.Reader(filename = vcffile)

    for chr in chrom_id_set:
        for record in vcf_reader.fetch(chr):

            chrom_id_set.add(record.CHROM)
            line_count['all'] = line_count.get('all', 0) + 1
            line_count[record.CHROM] = line_count.get(record.CHROM, 0) + 1

    return line_count, list(chrom_id_set)

def _set_step_num(line_count, sub_scale_num):

    step = line_count / sub_scale_num
    if step == 0: 
        step = 1
        sub_scale_num = line_count
    
    return sub_scale_num, step

if __name__ == '__main__':

    cmdopt = get_opt()
    main(cmdopt)
    print >> sys.stderr, '** For the flowers bloom in the desert **'
