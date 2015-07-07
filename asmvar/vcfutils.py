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
        
    def Add(self, mark, id, num, type, description):
        key = '##%s=<ID=%s' % (mark, id)
        val = ('##%s=<ID=%s,Number=%d,Type=%s,Description="%s">' % 
              (mark, id, num, type, description))
        self.header[key] = val
        return self

    def Record(self, headline):

        if   re.search (r'^##fileformat', headline): tag = '###'
        elif re.search (r'^#CHROM'      , headline): tag = '#CHROM'
        else: tag = headline.split(',')[0]
        self.header[tag] = headline

##
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

    def Add(self, key, context):
        self.info[key] = context
        return self
##
##
class Context(object): 

    def __init__(self):

        """
        VCF comtext
        """
        self.chrom = None
        self.pos   = None  
        self.Id    = None
        self.ref   = None  
        self.alt   = None 
        self.qual  = None    
        self.filters = None  
        self.info    = None   
        self.formats = None 
        self.sample  = []


        

