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
            # Stored all the index together as a array if the sequence are 
            # the same! But the hash key cannot in order as the sequence
            # The hash value is a list for recording the index of seq except
            # the tail(self.hashmer)
            self.hash_table[hash_id] = self.hash_table.get(hash_id, []) + [i]
			# This is a list to keep in order as sequence! And we should go 
            # thought the hash by scan this array to keep the right order!
            self.hash_pointer.append(hash_id)

