"""
===============================================
===============================================
Author: Shujia Huang
Date  : 2014-05-23 11:21:53
"""
import sys
import matplotlib.pyplot as plt

# My own class
import VariantDataManager as vdm
import VariantRecalibratorEngine as vre
import VariantRecalibratorArgumentCollection as VRAC

class VariantRecalibrator: 

    def __init__ (self):

        self.VRAC          = VRAC.VariantRecalibratorArgumentCollection()
        self.dataManager   = vdm.VariantDataManager()
        self.engine        = vre.VariantRecalibratorEngine(self.VRAC)
        self.badLodCutoff  = None
        self.LodCumInTrain = []
    def OnTraversalDone(self, data):

        self.dataManager.SetData(data)
        self.dataManager.NormalizeData()

        # Generate the positive model using the training data and evaluate 
        # each variant
        positiveTrainingData = self.dataManager.GetTrainingData()
        print >> sys.stderr, '[INFO] Training the goodModel ...'
        goodModel = self.engine.GenerateModel(positiveTrainingData, 
                                              self.VRAC.MAX_GAUSSIANS)
        print >> sys.stderr, ('[INFO] The converged information of goodModel '
                              'is:', goodModel.converged_)
        print >> sys.stderr, ('[INFO] The means of gaussion of goodModel is:\n', 
                              goodModel.means_)
        print >> sys.stderr, ('[INFO] The covariance of gaussion of goodModel '
                              'is:\n', goodModel.covars_, '\n')
        self.engine.EvaluateData(self.dataManager.data, goodModel, False)

        self.badLodCutoff, self.LodCumInTrain = self.dataManager.CalculateWorstLodCutoff()
        # Generate the negative model using the worst performing data and 
        # evaluate each variant contrastively
        print >> sys.stderr, '[INFO] Training the badModel ...'
        negativeTrainingData = self.dataManager.SelectWorstVariants(self.badLodCutoff)
        badModel = self.engine.GenerateModel(
            negativeTrainingData, 
            min(self.VRAC.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL, 
                self.VRAC.MAX_GAUSSIANS))

        print >> sys.stderr, ('\n[INFO] The converged information of badModel '
                              'is:' , badModel.converged_)
        print >> sys.stderr, ('[INFO] The means of gaussion of badModel '
                              'is:\n', badModel.means_)
        print >> sys.stderr, ('[INFO] The covariance of gaussion of badModel '
                              'is:\n', badModel.covars_, '\n')
        self.engine.EvaluateData(self.dataManager.data, badModel, True)

        if (not goodModel.converged_) or (not badModel.converged_): 
            raise ValueError ('[ERROR] NaN LOD value assigned. Clustering '
                              'with these variants and these annotations is '
                              'unsafe. Please consider raising the number of '
                              'variants used to train the negative model or '
                              'lowering the maximum number of Gaussians allowed '
                              'for use in the model.')

        # Find the VQSLOD cutoff values which correspond to the various 
        # tranches of calls requested by the user
        self.engine.CalculateWorstPerformingAnnotation(self.dataManager.data, 
                                                       goodModel, badModel)

    def VisualizationLodVStrainingSet(self, figName):

        fig = plt.figure()
        plt.title('LOD VS Positive training set', fontsize = 14)
        plt.plot(self.LodCumInTrain[:,0], self.LodCumInTrain[:,1], 'r-')
        #plt.scatter(self.LodCumInTrain[:,0], self.LodCumInTrain[:,1], c='r', 
        #            marker='.', linewidth = 0, alpha = 0.5)
        plt.plot([self.badLodCutoff, self.badLodCutoff], [0,1], 'g--')
        plt.ylim(0, 1.0)
        plt.xlim(-10, 10)
        plt.xlabel('Variant score threshold for the bad model', fontsize = 16)
        plt.ylabel('Rate of Positive->Negative', fontsize = 16)

        fig.savefig(figName + '.png')
        fig.savefig(figName + '.pdf')



