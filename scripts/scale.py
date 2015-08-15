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
import time
from subprocess import Popen, PIPE

import vcf
import pysam # I use the tabix module in pysam

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
    optp.add_option('--rec', dest = 'recursive', metavar = 'Boolen',
                    help = 'Recursvie to split the input vcffile(-v) by '
                    'its chromosome IDs(-c) or not.', default = '')

    # For Job running
    optp.add_option('-p', '--prefix', dest = 'outprefix', metavar = 'OUT',
                    help = 'Out file prefix', default = 'test')
    optp.add_option('--prog', dest = 'prog', metavar = 'STR',
                    help = 'The program for jobs', default = '')
    optp.add_option('--cmp', dest = 'comparameter', metavar = 'STR',
                    help = 'Common parameter for --prog', default = '')
    optp.add_option('--qsub', dest = 'qsub', metavar = 'QSUB',
                    help = 'qsub parameters for jobs', default = '')

    opt, _ = optp.parse_args()
    if not opt.vcffile: optp.error('Required [-v vcffile]\n')
    print >> sys.stderr, 'Parameters: python', ' '.join(sys.argv) 

    opt.number = abs(string.atoi(opt.number))
    opt.recursive = True if opt.recursive else False

    if opt.prog:
        opt.prog = os.path.abspath(opt.prog)

    # Create the outdir if it's not exists
    if not os.path.exists(opt.outdir):
        os.makedirs(opt.outdir)
    opt.outdir = os.path.abspath(opt.outdir)

    return opt

def main(opt):

    if opt.number > 0 or opt.recursive:
        sub_vcf_files = splitVCF(opt.vcffile, opt.ref_chrom, opt.number,
                                 opt.outdir + '/tmp_in_dir', opt.recursive)
    else:
        # Don't need to split the VCF file
        sub_vcf_files = [opt.vcffile]

    outinfo = createJobScript(opt.prog, opt.comparameter, sub_vcf_files, 
                              opt.outdir, opt.outprefix, opt.recursive)

    # Qsub the jobs
    if opt.qsub:
        qsubJobs(opt.qsub, [q[0] for q in outinfo])

    print '#Shell_Script\tOutput_file\tOutput_log'
    print '\n'.join(['\t'.join(s) for s in outinfo])

def qsubJobs(qsub_cmd, jobscripts):
    """
    Submitting jobs by qsub_cmd
    """
    import commands
    for q in jobscripts:

        sh_dir = os.path.dirname(os.path.abspath(q))
        (err_stat, job_id) = commands.getstatusoutput('cd %s && %s %s' % 
                                                      (sh_dir, qsub_cmd, q))
        if not err_stat:
            print >> sys.stderr, '[Good] Submitting job %s (%s) done' % (q, job_id)
        else:
            print >> sys.stderr, '[ERRR] Submitting job %s (%s) fail' % (q, job_id)
        

def createJobScript(program, com_parameters, input_files, 
                    outdir, outprefix, recursive):
    """
    Create job script.
    Args:
    """
    tmp_out_dir = outdir + '/tmp_out_dir'
    shell_dir   = outdir + '/shell'
    for d in (tmp_out_dir, shell_dir):
        if not os.path.exists(d):
            os.makedirs(d)

    outinfo   = []
    shell_pfx = shell_dir + '/' + outprefix
    for i, file in enumerate(input_files):

        fh = pysam.TabixFile(file)
        for chr in fh.contigs:

            sub_o_dir = tmp_out_dir + '/' + chr if recursive else tmp_out_dir
            if not os.path.exists(sub_o_dir):
                os.makedirs(sub_o_dir)
            outpfx = sub_o_dir + '/' + outprefix

            sub_out_file = '.'.join([outpfx, str(i + 1), chr, 'vcf'])
            sub_out_log  = '.'.join([outpfx, str(i + 1), chr, 'log'])
            sub_sh_file  = '.'.join([shell_pfx, str(i + 1), chr, 'sh'])
            outinfo.append([sub_sh_file, sub_out_file, sub_out_log])

            sh = open(sub_sh_file, 'w')
            sh.write('time python %s genotype %s -c %s -v %s > %s 2> %s\n' % 
                     (program, com_parameters, chr, file, sub_out_file,
                      sub_out_log))
            sh.close()

        fh.close()

    return outinfo

def splitVCF(vcffile, ref_chrom, split_num, sub_outdir, is_rec_split = True):

    print >> sys.stderr, '[INFO] ** Countting vcf lines. **'
    vcf_line_count, chrom_ids = _get_vcf_line_count(vcffile, ref_chrom)

    # get vcffile's file name by os.path.split
    _, fname = os.path.split(vcffile)
    if not os.path.exists(sub_outdir):
        os.makedirs(sub_outdir)

    print >> sys.stderr, '[INFO] ** Splitting vcf file. **'
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
            print >> sys.stderr, ('[INFO] ** Splitting VCF file of %s '
                                  'into %d. **' % (chrom, tot_num))
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

    print >> sys.stderr, '[INFO] ** Splited files all done. **'
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

        if line_num % 100000 == 0:
            print >> sys.stderr, ('[INFO] >> outputting %d lines in '
                                  'sub file. %s' % 
								  (line_num + 1, time.asctime()))

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

    print >> sys.stderr, ('[INFO] >> All outputting %d lines in '
                          'sub file. %s' % (line_num, time.asctime()))
    return sub_files

def _get_vcf_line_count(vcffile, chrom_id):
    """
    Get line number of this vcf file. And record them in a dict.
    """

    f = pysam.TabixFile(vcffile)
    chrom_id_set = set(f.contigs) # Initial
    f.close()

    if chrom_id:
        chrom_ids = set(chrom_id.split(','))

    line_count = {}
    vcf_reader = vcf.Reader(filename = vcffile)

    
    for chr in chrom_id_set:
        for record in vcf_reader.fetch(chr):

            line_count['all'] = line_count.get('all', 0) + 1
            line_count[record.CHROM] = line_count.get(record.CHROM, 0) + 1
            if line_count['all'] % 100000 == 0:
                print >> sys.stderr, ('[INFO] >> Countting %d lines. << %s' % 
                                      (line_count['all'], time.asctime()))

    print >> sys.stderr, '[INFO] ** The VCF line is %d' % line_count['all']
    return line_count, list(chrom_id_set)

def _set_step_num(line_count, sub_scale_num):

    step = line_count / sub_scale_num
    if step == 0: 
        print >> sys.stderr, ('[WARNING] The split file number is bigger than ' 
                              'the line number of VCF files. Reset it to be 1.')
        step = line_count
        sub_scale_num = 1

    if step * sub_scale_num < line_count:
        sub_scale_num += 1

    return sub_scale_num, step

if __name__ == '__main__':

    cmdopt = get_opt()
    main(cmdopt)
    print >> sys.stderr, '**>> For the flowers bloom in the desert <<**\n'
