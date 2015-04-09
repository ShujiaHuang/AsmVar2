"""
Common functions
"""
import os
import ctypes

dir = os.path.dirname(os.path.abspath(__file__))
encode = ctypes.CDLL(dir + '/encode.so')

class SeqHashTable(object):

    def __init__(self, seq, hashmer = None):
        """
        Initial hash
        """
        self.hashmer      = hashmer # The hash size, should be integer
        self.hash_table   = {}
        self.hash_pointer = []

        if hashmer is None: 
            return

        # Creat the sequence hash table and index
        self._set_hash(seq)
    
    def _set_hash(self, seq):
        """
        Build a hash table for the sequence.
        """

        if seq is None or len(seq) == 0 or len(seq) < self.hashmer:
            return;

        for i in range(len(seq) - self.hashmer + 1): 

            # Return an integer to represent the sequence
            hash_id = encode.hashEncode(seq[i:i + self.hashmer], self.hashmer)
            # Stored all the index together as an array if the sequence are 
            # the same! But the hash key cannot keep the order as the sequence.
            # Remember: The hash value is a list of the index of the seq except 
            # the tail(self.hashmer)
            self.hash_table[hash_id] = self.hash_table.get(hash_id, []) + [i]
			# This is a list to keep in order as sequence! And we should go 
            # thought the hash by scan this array to keep the right order!
            self.hash_pointer.append(hash_id)

def set_gap_open_penalty(seq, homopol_penalty):
    """
    Setting the gap open penalty by using the homopol_penalize model

    Args:
        ``length``: The haplotype sequence string.
        ``homopol_penalty``: penalty model could be ``CommonDatum.homopol_penalize``
                             and the ASCII value is from big to small.
    """

    gap_open_penalty = []
    # Fill in the penalty from the back(NOT from head!) to help 
    # left-justify indels
    homo_char = seq[-1].upper()
    hi = -1 # number of homopolymer 
    for c in seq[::-1]: # Count from tail

        if c.upper() == homo_char:
            if hi < len(homopol_penalty): # Not be overflow
                hi += 1 
        else:
            hi = 0

        penalty_value = ord(homopol_penalty[hi]) - 33
        if penalty_value < 0 or hi < 0:
            raise ValueError('Encountered negative gap open score: %s in %s '
                             '\n' % (penalty_value, seq))
        gap_open_penalty.append(chr(penalty_value)) # covert to be ASCII
        
        homo_char = c.upper() # Move homo_char to the next

    return ''.join(gap_open_penalty) # A char string











