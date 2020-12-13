"""
Simple classes designed to carry lots of attributes in a single instance.
"""
import re
import pprint
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
        self.acc = None
        self.name = None
        self.origin = None
        self.reverse = self.complement = self.revcomp = self.transcribe = None
        self.gff = self.protein = self.fasta = self.translate = self.mode = None
        self.gene = self.type = self.regexp = None
        self.match_field = self.match_value = None

        self.__dict__.update(kwds)

        # Parses out colon from data name if that exists.
        elems = self.acc.split(":")

        # Shortcut to types
        if len(elems) == 2:
            self.acc, self.type = elems

        # Shortcut to field/value matches.
        if len(elems) == 3:
            # Special casing the gene.
            self.acc, field, value = elems
            if field == "gene":
                self.gene = value
                self.type = "CDS"
            else:
                self.acc, self.match_field, self.match_value = elems

        self.start, self.end = utils.zero_based(start=self.start, end=self.end)
        self.regexp = re.compile(self.regexp) if self.regexp else None

    def unset(self):
        """
        Feature filtering parameters not set.
        """
        return not (self.start or self.end or self.type or self.gene or self.regexp)

    def pprint(self):
        pprint.pprint(self.__dict__)

    def __str__(self):
        return str(self.__dict__)