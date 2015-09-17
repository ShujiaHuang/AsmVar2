"""
=========================================
Variant quality score recalibratoe (VQSR)
=========================================
Author: Shujia Huang & Siyang Liu
Date  : 2014-05-23 11:21:53
"""
import sys
import re
import os
import sys
import time
import string
import optparse

# My own class
import VariantDataManager as vdm
import VariantRecalibrator as vror

def main(opt):

    # Just record the sites of training data 
    traningSet = vdm.LoadTrainingSiteFromVCF(opt.trainData)
    # Identify the traning sites
    hInfo, dataSet = vdm.LoadDataSet(opt.vcfInfile, traningSet)
    # init VariantRecalibrator object 
    vr = vror.VariantRecalibrator()
    # Traning modul and calculate the VQ for all the dataSet 
    vr.OnTraversalDone(dataSet)
    vr.VisualizationLodVStrainingSet(opt.figure + '.BadLodSelectInTraining')

    # For Record the Annnotations' values
    for d in vr.dataManager.annoTexts: 
        hInfo.add('INFO', d[0], 1, d[1], d[2])
    # Outputting the result as VCF format
    hInfo.add('INFO', 'VQ', 1, 'Float' , 'Variant Quality')
    hInfo.add('INFO', 'CU', 1, 'String', 'The annotation which was the worst '
              'performing in the Gaussian mixture modul, likely the reason why '
              'the variant was filtered out. It\'s the same tag as <culprit> '
              'in GATK')
    hInfo.add('INFO', 'NEGATIVE_TRAIN_SITE', 0, 'Flag', 'This variant was used '
              'to build the negative training set of bad variants')
    hInfo.add('INFO', 'POSITIVE_TRAIN_SITE', 0, 'Flag', 'This variant was used '
              'to build the positive training set of good variants')
    hInfo.add('INFO', 'InbCoeff', 1, 'Float' , 'Inbreeding coefficient '
              'as estimated from the genotype likelihoods per-sample when '
              'compared against the Hardy-Weinberg expectation')
 
    culprit, good, tot = {}, {}, 0.0
    annoTexts = [d[0] for d in vr.dataManager.annoTexts]
    idx = {c:i for i, c in enumerate(annoTexts)}

    for k, v in sorted(hInfo.header.items(), key = lambda d: d[0]): print v
    if opt.vcfInfile[-3:] == '.gz':
        I = os.popen('gzip -dc %s' % opt.vcfInfile)
    else:
        I = open(opt.vcfInfile)

    print >> sys.stderr, '\n[INFO] Outputting ...', time.asctime()
    n = 0
    j, monitor = 0, True
    while 1:

       lines = I.readlines(100000)
       if not lines: break

       for line in lines:

           n += 1
           if n % 100000 == 0: 
               print >> sys.stderr, '** Output lines %d %s' % (n, time.asctime())
           if re.search(r'^#', line): continue
           col = line.strip('\n').split()
           nratio = re.search(r';?NR=([^;]+)', col[7])
           if not nratio: continue

           atleastone = False
           for sample in col[9:]:
               field = sample.split(':')
               if field[0] == './.': continue
               atleastone = True
               break

           if not atleastone: continue
           
           order = col[0] + ':' + col[1]
           d  = dataSet[j]
           j += 1 # Increase the index of dataSet for the next cycle
           if d.variantOrder != order: 
               raise ValueError('[BUG] The order(%s) must be the same as '
                                'dataSet(%s)' % (order, d.variantOrder))

           # Deal with the INFO line
           vcfinfo = {}
           for info in col[7].split(';'): 
               k = info.split('=')[0]
               if monitor and k in vcfinfo: 
                   monitor = False
                   print >> sys.stderr, ('[WARNING] The tag: %s double hits in '
                                         'the INFO column at %s' % 
                                         (k, opt.vcfInfile))
               vcfinfo[k] = info

           tot += 1.0 # Record For summary
           culprit[annoTexts[d.worstAnnotation]] = culprit.get(
                annoTexts[d.worstAnnotation], 0.0) + 1.0 #For summary
           d.lod = round(d.lod, 2)
           for lod in [0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50]:
               if d.lod >= lod: good[lod] = good.get(lod, 0.0) + 1.0
        
           if d.atTrainingSite: 
                vcfinfo['POSITIVE_TRAIN_SITE'] = 'POSITIVE_TRAIN_SITE'
           if d.atAntiTrainingSite: 
                vcfinfo['NEGATIVE_TRAIN_SITE'] = 'NEGATIVE_TRAIN_SITE'

           vcfinfo['VQ'] = 'VQ=' + str(d.lod)
           vcfinfo['CU'] = 'CU=' + annoTexts[d.worstAnnotation]

           for text in annoTexts: 
               if text not in vcfinfo: 
                    vcfinfo[text] = (text + '=' + 
                                     str('%.2f' % d.raw_annotations[idx[text]]))

           col[7] = ';'.join(sorted(vcfinfo.values()))
           if d.lod < 0: d.lod = 0 # QUAL: donot allow value below 0
           col[5] = str(d.lod) # QUAL field should use phred scala

           print '\t'.join(col)
    I.close()
    print >> sys.stderr, '[INFO] Finish Outputting %d lines. %s' % (n, time.asctime())

    ## Output Summary
    print >> sys.stderr, '\n[Summmary] Here is the summary information:\n'
    for k, v in sorted(good.items(), key = lambda k:k[0]): 
        print >> sys.stderr, ('  ** Variant Site score >= %d: %d\t%0.2f' % 
            (k, v, v * 100 / tot))

    for k, v in sorted(culprit.items(), key = lambda k:k[0]):
        print >> sys.stderr, ('  ** Culprit by %s: %d\t%.2f' % 
            (k, v, v * 100.0 / tot))

def cmdopts():
    """
    The command line parameters for VQSR 
    """
    usage = ('\nUsage: %prog vqsr [--Train Training data set] '
             '[-i SampleVcfInfile] > Output')
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option('-i', '--InVcf', dest = 'vcfInfile', metavar = 'VCF', 
                    help = 'VCF for predict.', default = [])
    optp.add_option('-T', '--Train', dest = 'trainData', metavar = 'TRU', 
                    help = 'Traing data set at true  category', default = [])
    optp.add_option('-f', '--fig', dest = 'figure', metavar = 'FIG', 
                    help = 'The prefix of figure.', default = 'figtest')

    opt, _ = optp.parse_args()
    if len(opt.vcfInfile) == 0: optp.error('Required[-i vcfInfile]\n')
    if len(opt.trainData) == 0: optp.error('Required[-T trainData. VCFFormat]\n')
    print >> sys.stderr, ('[INFO] Parameters: python', sys.argv[0], 
                          '\n\t-i', opt.vcfInfile, 
                          '\n\t-T', opt.trainData, 
                          '\n\t-f', opt.figure, '\n')
    return opt

