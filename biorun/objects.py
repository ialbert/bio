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
        self.phase = "1"
        self.seqid = None
        self.gap_open = 11
        self.gap_extend = 1
        self.matrix = None
        self.json = None
        self.acc = None
        self.name = None
        self.inter = None
        self.reverse = self.complement = self.revcomp = self.transcribe = None
        self.gff = self.protein = self.fasta = self.translate = self.mode = None
        self.gene = self.type = self.regexp = None
        self.genome = False
        # Searching for attribute matches
        self.attr_name = self.attr_value = None

        self.__dict__.update(kwds)

        # Parses out colon from data name if that exists.
        elems = self.acc.split(":", 2)

        # Shortcut to types acc:type
        if len(elems) == 2:
            self.genome, self.fasta = False, True
            self.acc, self.type = elems

        # Shortcut to attrfield/value matches. acc:type:value
        if len(elems) == 3:
            self.genome, self.fasta = False, True
            self.acc, field, value = elems
            # Special casing the gene.
            if field == "gene":
                self.gene = value
                self.type = "CDS"
            else:
                self.acc, self.attr_name, self.attr_value = elems

        # An invalid parameter will be passed down as accession number.
        if self.acc.startswith("-"):
            msg = f"Unknown parameter: {self.acc}"
            utils.error(msg)

        # Allow commas in numbers, or sizes like 10Kb
        start = utils.parse_number(kwds.get("start", 0))
        end = utils.parse_number(kwds.get("end", 0))

        # Set zero based numbers.
        self.start, self.end = utils.zero_based(start=start, end=end)

        # Compile the regular expression.
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