"""
===============
Fusion vcf file
===============
Author: Shujia Huang
Date: 2015-09-29
"""
import os
import re
import sys
import optparse

dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir + '/../asmvar')
from utils import vcfutils

def main(opt):

    vcf_header = vcfutils.Header()
    get_new_header_record(opt.from_vcf, vcf_header, opt.add_format)
    get_new_header_record(opt.to_vcf, vcf_header) # Record header

    add_format_record, sam2idx = get_add_record(opt.from_vcf, 
                                                set(opt.add_format.split(':')))

    for k, h in sorted(vcf_header.header.items(), key = lambda d: d[0]): print h
    if opt.to_vcf[-3:] == '.gz': 
        I = os.popen('gzip -dc %s' % opt.to_vcf) 
    else: 
        I = open(opt.to_vcf)

    tot_num, use_num = 0, 0
    while 1:

        lines = I.readlines(100000)
        if not lines: break

        for line in lines:
            col = line.strip('\n').split()
            if re.search(r'^#CHROM', line):
                idx2sam = {i:sam for i, sam in enumerate(col[9:])}
            if line[0] == '#': continue

            tot_num += 1

            format_set = set(col[8].split(':')[1:])
            k = col[0] + ':' + col[1]
            if k not in add_format_record: continue

            use_num += 1

            fmat = {i:f for i, f in enumerate(col[8].split(':'))}
            for i, sam in enumerate(col[9:]):
                d = {fmat[i]:s for i, s in enumerate(sam.split(':'))}
                idx = sam2idx[idx2sam[i]]
                for f, v in add_format_record[k][idx].items():
                    format_set.add(f)
                    d[f] = str(v)

                format = [col[8].split(':')[0]] + sorted(list(format_set))
                col[i+9] = ':'.join(d[f] for f in format)

            col[8] = ':'.join(format)
            print '\t'.join(col)
    I.close()

    print >> sys.stderr, '[INFO] Total variants: %d; the number we left: %d' % (
        tot_num, use_num)

def get_add_record(vcffile, add_format_field_set):
    """
    """
    add_format_record = {} # 
    sam2idx = {}
    I = os.popen('gzip -dc %s' % vcffile) if vcffile[-3:] == '.gz' else open(vcffile)
    while 1:

        lines = I.readlines(100000)
        if not lines: break

        for line in lines:
            col = line.strip('\n').split()
            if re.search(r'^#CHROM', line):
                sam2idx = {sam:i for i, sam in enumerate(col[9:])}
            if line[0] == '#': continue

            fmat = {f:i for i, f in enumerate(col[8].split(':'))}
            k = col[0] + ':' + col[1]
            add_format_record[k] = []
            for i, sam in enumerate(col[9:]):
                format = sam.split(':')
                t = {}
                for f in (add_format_field_set & set(fmat.keys())):
                    t[f] = format[fmat[f]] if len(format) > fmat[f] else '.'
                add_format_record[k].append(t)

    I.close()

    return add_format_record, sam2idx
    
def get_new_header_record(vcffile, vcf_header, new_format_field = None):
    """
    Just read the header of vcffile and get the new record.
    Args:
        'vcffile': Input vcffile
        'new_format_field': 
    """
    new_format_fields = []
    if new_format_field:
        new_format_fields = ['##FORMAT=<ID=' + r for r in new_format_field.split(':')]

    I = os.popen('gzip -dc %s' % vcffile) if vcffile[-3:] == '.gz' else open(vcffile)
    for line in I:
        if line[0] != '#': break
        if new_format_fields:
            record = line.strip('\n').split(',')[0] # '##FORMAT<ID=%s'
            if record in new_format_fields:
                vcf_header.record(line.strip('\n'))
        else:
            # 'new_format_field' is empty, then record the header
            vcf_header.record(line.strip('\n'))
    I.close()
    

if __name__ == '__main__':

    usage = ('Usage: %prog --from [vcf] --add_format "<field1>:<field2>" '
             '--to [vcf] > output_vcf')
    optp = optparse.OptionParser(usage = usage)
    optp.add_option('--from', dest = 'from_vcf', metavar = 'VCF',
                    help = 'VCF for getting information.', default = '')
    optp.add_option('--to', dest = 'to_vcf', metavar = 'VCF',
                    help = 'VCF for adding information.', default = '')
    optp.add_option('--add_format', dest = 'add_format', metavar = 'STR',
                    help = 'Adding format field from "from_vcf" by parameter '
                    '"--from".', default = '')
    
    opt, _ = optp.parse_args()
    if len(opt.from_vcf) == 0: optp.error('Required[--from vcffile]\n')

    main(opt)
