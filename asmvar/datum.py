"""
This module contain all the usefull globle values for AsmVar2

The best thing is that ``datum`` should just be imported by ``executor``

"""
import os
import ctypes
import numpy as np

dir = os.path.dirname(os.path.abspath(__file__))
encode = ctypes.CDLL(dir + '/encode.so')

class Model(object):

    def __init__(self):
        """
        Initial
        """
        # Currently hard-coded insertion/deletion priors
        self.complex_deletion_prior = 5e-5
        self.complex_insertion_prior= 5e-6

        self.indel_prior_model = {1: "LIGC@:62/-*'&%$", 
                                  2: "LIGDB@><9630.,+**)(''&&%%%$$$", 
                                  3: "LIGA@B@><;8763220/.-,+++)*))(((''''&&&&&&%%%%%%%%$$$$$$$", 
                                  4: "LIGA@???=<886533210/.--,+**))))((('''''&&&&&&&&%%%%%%%%%%%$$$$$$$$", 
                                  5: 'LIGA@??>=>=;966543210///-,,++*', 
                                  6: 'LIGA@??>>=<=;:764532210/----,++', 
                                  7: 'LIGA@??>>==<;;987543210/....-,,,++++', 
                                  8: 'LIGA@??>>==<<;9876432200/..--,,,+++', 
                                  9: 'LIGA@??>>==<<;;9966432100//../..----,,,,,++++++', 
                                 10: 'LIGA@??>>==<<;;:986432110//..----,,,,++++', 
                                 11: 'LIGA@??>>==<<<;;:87642210////..--,,,,,+++', 
                                 12: 'LIGA@??>>==<<<;;;:986532110000/...-----,,,,,+++++', 
                                 13: 'LIGA@??>>==<<<;;;::987543111000/////.......--------,,,,,,,,,,,,,+++++++++', 
                                 14: 'LIGA@??>>==<<<;;;::987642210/0/.....-------,,,,,,,,+++++++', 
                                 15: 'LIGA@??>>==<<<;;;;::988754322110000////////.......------------,,,,,,,,,,,,,,,,,++++++++++', 
                                 16: 'LIGA@??>>==<<<;;;;:::98765321110////........-------,,,,,,,,,,,,,,+++++++++', 
                                 17: 'LIGA@??>>==<<<;;;;::::988764433211110000000///////.............-----------------,,,,,,,,,,,,,,,,,,,', 
                                 18: 'LIGA@??>>==<<<;;;:::::998875433221111000000///////.............-----------------,,,,,,,,,,,,,,,,,,,', 
                                 19: 'LIGA@??>>==<<<;;;;::::999887654433222221111111100000000//////////////..................------------', 
                                 20: 'LIGA@??>>==<<<;;;;::::9999876543322111000000///////............-----------------,,,,,,,,,,,,,,,,,,,', 
                                 21: 'LIGA@??>>==<<<;;;;::::9999988765544433322222221111111100000000000000//////////////////.............', 
                                 22: 'LIGA@??>>==<<<;;;;::::9999987765432221000000////////...........-----------------,,,,,,,,,,,,,,,,,,,', 
                                 23: 'LIGA@??>>==<<<;;;;::::9999998776543322111100000000////////................-------------------,,,,,,', 
                                 24: 'LIGA@??>>==<<<;;;;::::9999998887654433322111111100000000/////////////...................-----------'}



class CommonDatum(object):

    def __init__(self):

        ########## Common Value #############################################
        self.mot = -0.1 # use to onvert the map_qual to be a log10 value
        #####################################################################

        ########## Haplotype Alignment ######################################
        indel_base_error = [2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5,
                            1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3]
        extend_err = [1.4e-3 + 4.3e-4 * (n - 10) for n in range(11,50)]

        # errors score
        penalty = indel_base_error + extend_err
        # homopolymer indel error model and convert the value into ASCII
        self.homopol_penalty = ''.join([chr(int(33.5 + 10 * np.log((i + 1) * q) / np.log(0.1))) for i, q in enumerate(penalty)])

        # Alignment score matrix
        self.gap_extend = 3
        self.nucprior   = 2

        self.hashmer  = 7 # The unit size of hash
        self.hashsize = 4 ** self.hashmer
        self.hitN     = encode.hashEncode('N' * self.hashmer, self.hashmer)
        # The max allow size for aligniing process. 
        # Set this value just want to limit the time
        #self.max_align_size = self.hashsize 
        self.max_align_size = 2000 

        self.use_read_mapq = True # Decide by myself and used in `alignment`
        self.do_calcu_flank_score = False # Decide by myself used in `alignment`
        #####################################################################

        ########## Some common values #######################################
        self.min_float  = 1e-300 # The min float value.
        ########## Genotyping ###############################################
        ## For E-M step
        self.max_iter_num = 50
        #####################################################################


