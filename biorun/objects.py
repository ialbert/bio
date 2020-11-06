"""
Simple classes designed to carry lots of attributes in a single instance.
"""
import re
from biorun import utils, const


class Param(object):
    """
    Stores all the input parameters that have grown too numerous to pass individually.
    """

    def __init__(self, **kwds):
        self.start = self.end = 0
        self.phase = "1"
        self.seqid = None
        self.gap_open = 11
        self.gap_extend = 1
        self.matrix = None
        self.json = None
        self.name = None
        self.reverse = self.complement = self.revcomp = self.transcribe = None
        self.gff = self.protein = self.fasta = self.translate = self.mode = None
        self.gene = self.type = self.regexp = None

        self.__dict__.update(kwds)

        # Parses out colon from data name if that exists.
        if self.name and ":" in self.name:
            self.name, word = self.name.split(":")
            self.gene, self.type = word, 'CDS'

        self.start, self.end = utils.zero_based(start=self.start, end=self.end)
        self.regexp = re.compile(self.regexp) if self.regexp else None

    def unset(self):
        """
        Feature filtering parameters not set.
        """
        return not (self.start or self.end or self.type or self.gene or self.regexp)

    def __str__(self):
        return str(self.__dict__)