"""
This module contain fucntions and class for genotyping the variats.

This module will connect the `Haplotype` and `alignment` class
"""
import copy
from itertools import combinations as itertools_combinations
import numpy as np

import time


import vcf
import alignment as alg # The alignment module
from haplotype import Haplotype

import datum as DM
COMDM = DM.CommonDatum()

class Diploid(object):
    """
    A class represent a diploid genotype. It store two haplotype
    """
    def __init__(self, haplotype1, haplotype2):
        """
        Constructor.

        Args:
            `haplotype1`: The data type is Haplotype
            `haplotype2`: The data type is Haplotype
        """
        if haplotype1.hapregion_hash_id != haplotype2.hapregion_hash_id:
            raise ValueError('[ERROR] The two haplotypes did not have the same '
                             'haplotype region, could not build up a diploid.')

        self.hap1 = haplotype1 # Do I have use deepcopy here?
        self.hap2 = haplotype2 # Do I have use deepcopy here?

    def _cal_realign_regions(self):
        """
        """
        boundary = 1
        regions  = set()
        if len(self.hap1.variants) == 0 and len(self.hap1.variants) == 0:
            # Small haplotype or it may be a big haplotype but it's a refernce
            # haplotype and have no variants in it and just cut the load reads
            # region to fix the max aligment size.
            start = self.hap1.hap_start - 1
            end   = min(self.hap1.hap_start + COMDM.max_align_size, 
                        self.hap1.hap_end)
            regions.add((start, end))
        else: 
            # Hete or homo variants genotype.
            for v in self.hap1.variants:
                regions.add((max(0, v.POS - boundary), v.POS + boundary))
            for v in self.hap2.variants:
                regions.add((max(0, v.POS - boundary), v.POS + boundary))
            """
            if len(self.hap1) < COMDM.max_align_size:
                regions.add((self.hap1.hap_start - 1, self.hap1.hap_end))
            else:
                # Big haplotype! We should not load all the reads in this 
                # process or it'll cost a lot of time in mapping. We just 
                # prefer reads which loaded around the variants' breakpoints. 
                for v in self.hap1.variants:
                    regions.add((max(0, v.POS - boundary), v.POS + boundary))

            if len(self.hap2) < 0:
                regions.add((self.hap2.hap_start - 1, self.hap2.hap_end))
            else:
                for v in self.hap2.variants:
                    regions.add((max(0, v.POS - boundary), v.POS + boundary))
            """

        # Cause: Some regions in 'regions' may overlap with others
        # but I'll deal with this situation in `alg.alignReadToHaplotype`
        return self._merge(sorted(list(regions)))

    def _merge(self, regions, dis_delta = 1):
        """
        Merge the region if them happen to overlap with others
        retrun a list of sorted regions.

        Args:
            'regions': A list of sorted regions. Format: (start, end)
        """
        new_reg = []
        pre_pos = regions[0][0] # pre start postion
        flag    = False
        reg_start, reg_end = 0, 0
        for start, end in regions:
            if end < start:
                raise ValueError('The region start > end (%d, %d). This is '
                                 'not allow when call merge function in ge-'
                                 'notype.' % (start, end))
            if not flag: # The light is on => get region
                reg_start, reg_end = start, end
                flag = True
            else:
                if pre_pos > start:
                    raise ValueError("Your regions haven't been sorted.\n")
                if reg_end + dis_delta >= start:
                    # Overlap and should be merged
                    if end > reg_end: reg_end = end
                else:
                    new_reg.append((reg_start, reg_end))
                    reg_start, reg_end = start, end
            pre_pos = start

        if flag: new_reg.append((reg_start, reg_end))

        return new_reg
        
    def calLikelihood(self, read_buffer_dict, bam_reader):
        """
        Calculate the genotype likelihood for this single bam_reader's sample.

        Args:
            `read_buffer_dict`: A hash to reads. It's a contianer of reads, 
                                use it to prevent recalculting the hash-sequence 
                                for the same reads. This could save the running 
                                time.
            `bam_readers`: A single bamfile reader opened by `pysam.AlignmentFile`

        return a log10 likelihood.
        """
        # List of likelihood for each aligning read. each value in the array
        # represent a likelihood value of read align to the haplotype
        # Remember these likelihood will be changed follew different bamfile
        regions = self._cal_realign_regions()
        self.hap1.likelihood = alg.alignReadToHaplotype(self.hap1,
                                                        read_buffer_dict,
                                                        bam_reader,
                                                        regions)
        self.hap2.likelihood = alg.alignReadToHaplotype(self.hap2,
                                                        read_buffer_dict,
                                                        bam_reader,
                                                        regions)
        lksize1 = len(self.hap1.likelihood)
        lksize2 = len(self.hap2.likelihood)
        if lksize1 != lksize2:
            raise ValueError('[ERROR] The two haplotype in this diploid must'
                             'have the same number alignment of reads. But '
                             'the number is %s and %s .' % (lksize1, lksize2))

        likelihood = None
        if lksize1: # Not empty 
            likelihood = 0.0

        # Calculate the log10 likelihood for this diploid region
        for i in range(lksize1):

            log10lk1 = self.hap1.likelihood[i] 
            log10lk2 = self.hap2.likelihood[i] 

            # Almost good to 1000 times. Just take the highest and forget 
            # the small one
            if abs(log10lk1 - log10lk2) >= 3.0:
                likelihood += (np.log10(0.5) + max(log10lk1, log10lk2))

            # The likelihoods of the two haplotypes are almost the same. This
            # could just happens when either they both bad fits or both are 
            # perfect match. Any way the read lies in a position that it 
            # cannot be used to distingush the two haplotypes.
            elif abs(log10lk1 - log10lk2) <= 1e-3:
                likelihood += log10lk1

            # Calculate as the normal way: Combine the two likelihood
            else:
                prob = 0.5 * (np.power(10, log10lk1) + np.power(10, log10lk2))
                likelihood += np.log10(prob) # now it's a log10 value

        # Generally, for our situation the likelihood should always < 0, but 
        # it may > 0, once we adjust with flank region of haplotype in the 
        # `alignReadToHaplotype` process.
        read_map_count = lksize1 # Count of mapping reads of this sample

        best_loglikelihood = -1e20
        return max(best_loglikelihood, likelihood), read_map_count


def generateAllGenotypes(haplotypes):
    """
    Generate a list of potentail genotypes. 
    """
    gentypes   = []
    index = range(len(haplotypes))
    for i in index:
        for j in index[i:]:
            # Create Diploid
            # Each element is [diploid, hap_hash_id1, hap_hash_id2]
            # hash_id is the identify code of each haplotype!
            gentypes.append([Diploid(haplotypes[i], haplotypes[j]), 
                             hash(haplotypes[i]), hash(haplotypes[j])])

    return gentypes

def generateAllHaplotypeByVariants(ref_fa_stream, max_read_len, winvar):
    """
    Generate all the potential haplotype in the window of `winwar`.
    And return a list of haplotype, corresponing to all the possible
    combination of variants in this specific window

    Args:
        `ref_fa_stream`: Input fasta stream of reference
        `max_read_len`: The max length of reads
        `winvar`: It's dict(chrom=chrom, start=start, end=end, variant=var)
    """
    # This is the reference haplotype attempt that no variants in this window
    #refhap = Haplotype(ref_fa_stream, winvar['chrom'], winvar['start'], 
    #                   winvar['end'], max_read_len)
    haplist = [] # Initial with empty list
    num_var = len(winvar['variant'])
    vindex  = range(num_var) # An array contain all the index of variants
    var_overlap_index = _get_overlap_var_index(winvar['variant'])
    # Generate all the combination of haplotypes
    # All the combinantion of haplotype may be too much! 

    done = set() # A set to save the done variants' combination
    for n in range(num_var):
        n += 1
        for idxs in itertools_combinations(vindex, n):

            new_idxs = _recreate_varaint_idxlist_by_overlop(idxs,
                                                            winvar['variant'],
                                                            var_overlap_index)
            for idx in new_idxs:

                hash_id = hash(tuple(idx))
                if hash_id in done: continue # Ignore all the done part!
                done.add(hash_id)

                # Must be deep copy, prevent to changing the raw `winvar`
                varlist = [copy.deepcopy(winvar['variant'][j]) for j in idx]
                need_iter = True
                while need_iter:
                    """
                    Iter all the possible variant combination, and each
                    variant combination could be join to be one haplotype
                    """
                    need_iter = False
                    var = []
                    for i, v in enumerate(varlist):

                        tmp = copy.deepcopy(v) # Still have to deep copy
                        if len(v.ALT) > 1:
                            # v.ALT array size will be smaller here
                            need_iter = True #  
                            tmp.ALT   = [v.ALT.pop()]
                        var.append(copy.deepcopy(tmp))

                    hap = Haplotype(ref_fa_stream, winvar['chrom'],
                                    winvar['start'], winvar['end'], 
                                    max_read_len, var)
                    haplist.append(hap)

    return haplist

def _get_overlap_var_index(varants):
    """
    Self overlap. Find all the overlap relationship with other element
    in 'variants' and return a hash which record the overlap index for
    all the array's elements
    """
    var_overlap_index = {}
    all_indexs = range(len(varants))
    loop_i     = 1
    for idx in all_indexs:

        start = varants[idx].POS
        end   = varants[idx].POS + len(varants[idx].REF) - 1
        flag  = True
        var_overlap_index[idx] = [idx] # Initial with the own index

        for i in all_indexs[loop_i:]:
            if idx == i: continue
            if start > varants[i].POS + len(varants[i].REF) - 1: continue
            if end < varants[i].POS: break

            # Overlap now!
            if flag:
                loop_i = i
                flag   = False

            var_overlap_index[idx].append(i)

    return var_overlap_index

def _recreate_varaint_idxlist_by_overlop(index_list, 
                                         variants,
                                         index_overlap_hash): 
    """
    Separating the overlap variants' index into different haplotype.
    
    return the new index list and it's a 2-D array.
    ####
    Args:
        `index_list`: A list of array's indexs
        `index_overlap_hash`: A hash contain the overlap relationship
    """
    pre_index   = index_list[0]
    index_nodes = [[pre_index]] # The first index
    for i in index_list[1:]:
        # Getting Overlap and creat index nodes
        if i in index_overlap_hash[pre_index]:
            # Collecte all the overlap variant's index into one node 
            index_nodes[-1].append(i)
        else:
            index_nodes.append([i]) # New node

        if (variants[pre_index].POS + len(variants[pre_index].REF) - 1 < 
            variants[i].POS + len(variants[i].REF) - 1):
            pre_index = i

    new_index_list = []
    need_iter = True
    while need_iter:

        need_iter = False
        tmp_index_list = []
        for index in index_nodes:
            # get index by looping all the index nodes.
            i = index[0] # default
            if len(index) > 1:
                need_iter = True
                i = index.pop()

            tmp_index_list.append(i)

        new_index_list.append(tmp_index_list)

    return new_index_list

