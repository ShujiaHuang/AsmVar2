"""
================================================
My own Gaussion Mixture Model for SV genotyping.
Learn form scikit-learn
================================================

Author: Shujia Huang
Date  : 2014-01-06 14:33:45

"""
import sys
import string
import re
import os
import gzip

import numpy as np
import scipy.stats as sp_stats
from sklearn.metrics import roc_curve

# My own class
dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir + '/..')
from utils import vcfutils
import VariantRecalibratorArgumentCollection as VRAC
import VariantDatum as vd

class VariantDataManager:

    def __init__(self, data=None):
        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        self.annotationMean = None
        self.annotationSTD  = None
        self.annoTexts      = [['InbCoeff', 'Float', 'Inbreeding coeff-'
                                'icient as estimated from the genotype likeli-'
                                'hoods per-sample when compared against the '
                                'Hardy-Weinberg expectation'], \
                               ['NRatio', 'Float', 'N ratio of ALT sequence'],\
                               ['HR', 'Integer', 'Homozygous run'], \
                               ['FS', 'Float', 'Phred-scaled p-value using '
                                'Fisher\'s exact test to detect strand bias']]

        self.data = [] # list < VariantDatum >
        if data: # data is not None
            if not isinstance(data[0],vd.VariantDatum): 
                raise ValueError('[ERROR] The data type should be '
                                 '"VariantDatum" in VariantDataMa-'
                                 'nager(),but found %s'% str(type(data[0])))
            self.data = data
            for i, d in enumerate(self.data): 
                self.data[i].annotations = np.array(self.data[i].annotations)

    def SetData(self, data):

        if not isinstance(data[0],vd.VariantDatum): 
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in VariantDataManager(),but found %s' % 
                             str(type(data[0])))
        self.data = data
        for i, d in enumerate(self.data): 
            self.data[i].annotations = np.array(d.annotations)

    def NormalizeData(self):

        data = np.array([d.annotations for d in self.data], dtype=float)
        mean = data.mean(axis=0); self.annotationMean = mean
        std  = data.std(axis=0) ; self.annotationSTD  = std

        # foundZeroVarianceAnnotation
        if any(std < 1e-5): 
            raise ValueError('[ERROR] Found annotations with zero variance. '
                             'They must be excluded before proceeding.')
        
        # Each data now is (x - mean)/sd
        for i, d in enumerate(data):
            self.data[i].annotations = (d - mean) / std 
            # trim data by standard deviation threshold and mark failing data 
            # for exclusion later
            self.data[i].failingSTDThreshold = False
            if any(np.abs(self.data[i].annotations) > self.VRAC.STD_THRESHOLD):
                self.data[i].failingSTDThreshold = True

    def GetTrainingData(self):

        trainingData = [d for d in self.data if ((not d.failingSTDThreshold) 
                                                 and d.atTrainingSite)]
        print >> sys.stderr, ('[INFO] Training with %d variants after standard '
                              'deviation thresholding.' % len(trainingData))

        if len(trainingData) < self.VRAC.MIN_NUM_BAD_VARIANTS:
            print >> sys.stderr, ('[WARNING] Training with very few variant '
                                  'sites! Please check the model reporting '
                                  'PDF to ensure the quality of the model is '
                                  'reliable.')
        if len(trainingData) > self.VRAC.MAX_NUM_TRAINING_DATA:
            print >> sys.stderr, ('[WARING] Very large training set detected. '
                                  'Downsampling to %d training variants.' % 
                                  self.VRAC.MAX_NUM_TRAINING_DATA)
            np.random.shuffle(trainingData) # Random shuffling
            return list(trainingData[i] 
                        for i in range(self.VRAC.MAX_NUM_TRAINING_DATA))

        return trainingData 

    def SelectWorstVariants(self, badLod):

        trainingData = []
        for i,d in enumerate(self.data):
            if(d.lod < badLod) and (not d.failingSTDThreshold):
                trainingData.append(d)
                self.data[i].atAntiTrainingSite = True # I do need: i order to be the same as self.data 
        print >> sys.stderr, '[INFO] Training with worst %d scoring variants --> variants with LOD < %.2f.' %(len(trainingData), badLod)

        if len(trainingData) > self.VRAC.MAX_NUM_TRAINING_DATA:
            print >> sys.stderr, '[WARING] Very large training set detected. Downsampling to %d training variants.' % self.VRAC.MAX_NUM_TRAINING_DATA
            np.random.shuffle(trainingData) # Random shuffling
            return list(trainingData[i] for i in range(self.VRAC.MAX_NUM_TRAINING_DATA))

        return trainingData

    def CalculateWorstLodCutoff(self):

        lodThreshold, lodCum = None, []
        if len(self.data) > 0:

            lodDist = np.array([[d.atTrainingSite, d.lod] 
                                for d in self.data if(not d.failingSTDThreshold)])

            # I just use the 'roc_curve' function to calculate the worst 
            # LOD threshold, not use it to draw ROC curve And 'roc_curve' 
            # function will output the increse order, so that I don't 
            # have to sort it again
            _, tpr, thresholds = roc_curve(lodDist[:,0], lodDist[:,1]) 
            lodCum = [[thresholds[i], 1.0 - r] for i, r in enumerate(tpr)]

            for i, r in enumerate(tpr):
                if r > 1.0 - self.VRAC.POSITIVE_TO_NEGATIVE_RATE: 
                    lodThreshold = round(thresholds[i])
                    break

        return lodThreshold, np.array(lodCum)

def SampleFaLen(faLenFile):

    if faLenFile[-3:] == '.gz': 
        I = os.popen('gzip -dc %s' % faLenFile)
    else: 
        I = open(faLenFile)
    
    data = {}
    while 1:
        lines = I.readlines(100000)
        if not lines: break
        for line in lines:
            col = line.strip('\n').split()
            data[col[0]] = string.atoi(col[1])
    I.close()
    return data

def LoadTrainingSiteFromVCF(vcffile):
    """
    Just record the training site positions
    """

    if vcffile[-3:] == '.gz':
        I = os.popen('gzip -dc %s' % vcffile)
    else:
        I = open(vcffile)

    dataSet = set()
    while 1:

        lines = I.readlines(100000)
        if not lines: break

        for line in lines:

            if re.search(r'^#', line): continue
            col = line.strip('\n').split()
            dataSet.add(col[0] + ':' + col[1])
    I.close()

    return dataSet

def LoadDataSet(vcfInfile, traningSet):

    if len(traningSet) == 0: 
        raise ValueError('[ERROR] No Training Data found')

    I = gzip.open(vcfInfile) if vcfInfile[-3:] == '.gz' else open(vcfInfile)
    data, hInfo = [], vcfutils.Header()
    while 1: # VCF format

        lines = I.readlines(100000)
        if not lines: break

        for line in lines:

            col = line.strip('\n').split()
            if re.search(r'^#CHROM', line): 
                col2sam = {i+9:sam for i,sam in enumerate(col[9:])}

            # Record the header information
            if re.search(r'^#', line):
                hInfo.record(line.strip('\n'))
                continue

            # Get inbreeding coefficient. If fail then recalculating.
            # It's calculated like: 1.0 - hetCount/Expected_hetCount in VCF
            inbCoeff = re.search(r';?InbCoeff=([^;]+)', col[7])
            if not inbCoeff:
                # We have to calculate the ibreedCoeff if there's no one
                inbCoeff = vcfutils.calcuInbreedCoeff([c.split(':')[0] 
                                                       for c in col[9:]])
            else:
                inbCoeff = float(inbCoeff.group(1))
            inbCoeff = round(inbCoeff, 2)

            nratio = re.search(r';?NR=([^;]+)', col[7])
            nratio = round(float(nratio.group(1)), 2)

            hom_run = re.search(r';?HR=([^;]+)', col[7])
            hom_run = round(float(hom_run.group(1)), 2)

            fmat = {k:i for i, k in enumerate(col[8].split(':'))} # Get Format
            for tag in ['SB']:
                if tag not in fmat: raise ValueError('[ERROR] The "Format" '
                                                     'fields did not contian '
                                                     '"%s" in VCF: %s\nAT: %s\n' 
                                                     %(tag, vcfInfile, line))

            fs = 0 # The phred scale of strand bias's pvalue
            atleastone = False
            for sample in col[9:]: 

                field = sample.split(':')
                if field[0] == './.': continue
                atleastone = True

                sb = [int(d) for d in field[fmat['SB']].split(',')]
                if len(sb) != 4: 
                    raise ValueError('[ERROR] "SB" must just get 4 '
                                     'elements.\n')
                ##################################################
                # Fisher exact test:                             #
                #              [ref-fwd,     ref-reverse],       #
                #              [non-ref-fwd, non-ref-reverse]    #
                ##################################################
                _, pvalue = sp_stats.fisher_exact([[sb[0], sb[1]], 
                                                   [sb[2], sb[3]]])
                fs += -10 * np.log10(pvalue) # Sum up

            if not atleastone: continue 

            datum = vd.VariantDatum()
            datum.raw_annotations = [inbCoeff, nratio, hom_run, fs]
            datum.annotations  = [inbCoeff, nratio, hom_run, fs]
            datum.variantOrder = col[0] + ':' + col[1]
            if datum.variantOrder in traningSet: 
                datum.atTrainingSite = True

            data.append(datum)

    I.close()

    return hInfo, np.array(data)

