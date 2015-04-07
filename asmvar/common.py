"""
Common functions
"""
import os
import ctypes

dir = os.path.dirname(os.path.abspath(__file__))
encode = ctypes.CDLL(dir + '/encode.so')

class SeqHashTable(object):

    def __init__(self, seq, hashmer):
        """
        Initial hash
        """
        self.hashmer      = hashmer # The hash size, should be integer
        self.hash_table   = {}
        self.hash_pointer = []

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
            #print '** i =', i, 'hash_id =', hash_id, 'hashmer =', self.hashmer, 'seq =', seq[i:i + self.hashmer]
            # Stored all the index together as a array if the sequence are 
            # the same! But the hash key cannot in order as the sequence
            self.hash_table[hash_id] = self.hash_table.get(hash_id, []) + [i]
            self.hash_pointer.append(hash_id) # This is in order as sequence!

