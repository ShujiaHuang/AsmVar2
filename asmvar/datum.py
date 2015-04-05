"""
This module contain all the usefull globle values for AsmVar2
"""

class Genotype(object):

    def __init__(self):
        base_error = [2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5,
                      1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3]
        extend_err = [1.4e-3 + 4.3e-4 * (n - 10) for n in range(11,50)]

        # Indel errors for each base
        self.base_indel_error = base_error + extend_err
        # Convert the base indel error value into ASCII
        self.indel_error_qual = ''.join([chr(int(33.5 + 10 * log((i + 1) * q) / log(0.1)))
                                for i, q in enumerate(self.base_indel_error)])

 class Align(object):

    def __init__(self):
        self.hashmer  = 7 # The unit size of hash
        self.hashsize = 4 ** self.hashmer
        # The max allow size for aligniing process. 
        # Set this value just want to limit the time
        self.max_align_size = self.hashsize 

