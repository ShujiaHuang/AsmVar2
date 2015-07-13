"""
A class for output VCF file. PyVCF does not able to add or update information 
fields for sample's FORMAT field. That make us have to create another classes
(like these) to handle that problem
"""
import re

class Header(object):

    def __init__(self, hInfo = None): 
        """
        VCF header information
        """
        self.header = {}
        if hInfo and (type(hInfo) is not dict): 
            raise ValueError ('The data type should be "dict" in class '
                              'of "VCFHeader", but found %s' % str(type(hInfo)))
        if hInfo: self.header = hInfo
        
    def add(self, mark, id, num, type, description):
        key = '##%s=<ID=%s' % (mark, id)
        val = ('##%s=<ID=%s,Number=%s,Type=%s,Description="%s">' % 
              (mark, id, num, type, description))
        self.header[key] = val
        return self

    def record(self, headline):

        if   re.search (r'^##fileformat', headline): tag = '###'
        elif re.search (r'^#CHROM'      , headline): tag = '#CHROM'
        else: tag = headline.split(',')[0]
        self.header[tag] = headline


class Info(object): 

    def __init__(self, info = None):
        """
        INOF fields information
        """
        self.info = {}
        if info and (type(info) is not dict): 
            raise ValueError ('The data type should be "dict" in class '
                              'of "VCFInfo", but found %s' % str(type(info)))
        if info: self.info = info

    def add(self, key, context):
        self.info[key] = context
        return self


class Context(object): 

    def __init__(self):

        """
        VCF comtext
        """
        self.chrom = None
        self.pos   = None  
        self.Id    = None
        self.ref   = None  
        self.alt   = []
        self.qual  = None    
        self.filter = []
        self.info   = {} 
        self.format = None 
        self.sample = []

    def print_context(self):
        """
        """
        if self.chrom:

            print '\t'.join([self.chrom,
                             str(self.pos),
                             '.' if not self.Id else self.Id,
                             self.ref,
                             ','.join(self.alt),
                             str(self.qual),
                             '.' if not self.filter else ','.join(self.filter),
                             '.' if not self.info else ';'.join(','.join(v) 
                                             for k, v in self.info.items()),
                             ':'.join(self.format),
                             '\t'.join(self.sample)])



        

