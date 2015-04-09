"""
This module contain all the usefull globle values for AsmVar2

The best thing is that ``datum`` should just be imported by ``executor``

"""
import numpy as np

class CommonDatum(object):

    def __init__(self):
        indel_base_error = [2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5,
                            1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3]
        extend_err = [1.4e-3 + 4.3e-4 * (n - 10) for n in range(11,50)]

        # errors score
        penalty = indel_base_error + extend_err
        # homopolymer indel error model and convert the value into ASCII
        self.homopol_penalty = ''.join([chr(int(33.5 + 10 * np.log((i + 1) * q) / np.log(0.1)))
                                        for i, q in enumerate(penalty)])

        # Alignment score matrix
        self.gap_extend = 3
        self.nucprior   = 2

        self.hashmer  = 7 # The unit size of hash
        self.hashsize = 4 ** self.hashmer
        # The max allow size for aligniing process. 
        # Set this value just want to limit the time
        self.max_align_size = self.hashsize 

