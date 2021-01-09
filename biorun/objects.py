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

        # List all values to help editor figure out what are valid.
        self.acc = ''
        self.db = ''
        self.phase = "1"
        self.seqid = None
        self.gap_open = 11
        self.gap_extend = 1
        self.matrix = None
        self.json = None
        self.uid = None
        self.name = None
        self.inter = None
        self.update = False
        self.reverse = self.complement = self.revcomp = self.transcribe = None
        self.gff = self.protein = self.fasta = self.translate = self.mode = None
        self.gene = self.type = self.regexp = None

        # Output data should be a record type (as opposed to a sequence).
        self.record = False

        # Output features rather than genome.
        self.features = False

        # Searching for attribute matches
        self.attr_name = self.attr_value = None

        # Override values from the constructor.
        self.__dict__.update(kwds)

        # Parses out colon from data name if that exists.
        elems = self.acc.split(":", 1)

        # The first element is the accession.
        self.acc = elems[0]

        # Shortcuts: acc:type or acc:key=value
        if len(elems) == 2:

            # We are selecting for features.
            self.features = True

            # Parse out the components
            tmp, text = elems

            # See if the second filed is splittable
            pieces = text.split("=", 1)

            # The "shorter" shortcut.
            if len(pieces) == 1:
                self.gene = text
                self.type = "CDS"
            else:
                field, value = pieces
                # A few fields are special cased.
                if field.lower() == "name":
                    self.name = value
                elif field.lower() == "id":
                    self.uid = value
                else:
                    self.attr_name, self.attr_value = field, value

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
