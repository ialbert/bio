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
        self.acc = ''
        self.phase = "1"
        self.seqid = None
        self.gap_open = 11
        self.gap_extend = 1
        self.matrix = None
        self.json = None
        self.uid = None
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
        elems = self.acc.split(":", 1)

        # The first element is the accession.
        self.acc = elems[0]

        # Shortcut to types acc:type
        if len(elems) == 2:
            tmp, text = elems

            # See if the second filed is splittable
            pieces = text.split("=", 1)

            # Shortcuts:
            # acc:foo ==> --name foo, --type CDS
            # acc:key=value ==> matches attribute
            if len(pieces) == 1:
                self.gene = text
            else:
                field, value = pieces
                # A few fields are not attributes thus handled differently.
                if field.lower() == "name":
                    self.name = value
                elif field.lower() == "id":
                    self.uid = value
                else:
                    self.attr_name, self.attr_value = field, value

        # When setting genes select for coding sequences.
        if self.gene:
            self.type = "CDS"

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