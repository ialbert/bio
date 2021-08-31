import re
import sys

try:
    from Bio import Entrez
except ImportError as exc:
    print(f"# Error: {exc}", file=sys.stderr)
    print(f"# This program requires biopython", file=sys.stderr)
    print(f"# Install: conda install -y biopython>=1.79", file=sys.stderr)
    sys.exit(-1)

from biorun.libs import placlib as plac
from tqdm import tqdm
from biorun import utils
from urllib.error import HTTPError

Entrez.email = 'foo@foo.com'

ncbi_patt = r'(?P<letters>[a-zA-Z]+)(?P<under>_?)(?P<digits>\d+)(\.(?P<version>\d+))?'
ncbi_patt = re.compile(ncbi_patt)


#
# Refseq accession numbers
#
# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
#

class Accession:
    """
    Represents an accession number
    """

    def __init__(self, value):
        self.value = value

        try:
            int(self.value)
            self.is_int = True
        except ValueError as exc:
            self.is_int = False

        m = ncbi_patt.search(value)

        # Valid accession numbers
        self.valid = bool(m) or self.is_int

        self.code = m.group("letters") if m else ''
        self.digits = m.group("digits") if m else ''
        self.under = m.group("under") if m else ''
        self.version = m.group("version") if m else ''

        # It is refseg if it has underscore
        self.refseq = bool(self.under)

    def is_nucleotide(self):
        if self.code:
            cond = self.code in ["AC", "NC", "NG", "NT", "NW", "NZ", "NM", "XM", "XR"]
        else:
            num1, num2 = len(self.code), len(self.digits)
            cond = (num1 == 1 and num2 == 5) or (num1 == 2 and num2 == 6) or (num1 == 3 and num2 == 8)

        return cond

    def is_protein(self):
        if self.refseq:
            cond = self.code in ["AP", "NP", "YP", "XP", "WP"]
        else:
            num1, num2 = len(self.code), len(self.digits)
            cond = (num1 == 3 and num2 == 5) or (num1 == 3 and num2 == 7)
        return cond


def efetch(ids, db, rettype='gbwithparts', retmode='text'):
    try:
        stream = Entrez.efetch(db=db, id=ids, rettype=rettype, retmode=retmode)
        stream = tqdm(stream, unit='B', unit_divisor=1024, desc='# downloaded', unit_scale=True, delay=5, leave=False)
    except HTTPError as exc:
        utils.error(f"Accession or database may be incorrect: {exc}")

    for line in stream:
        print(line, end='')
        stream.update(len(line))
    stream.close()


@plac.pos("acc", "accession numbers")
@plac.opt("db", "database", choices=["nuccore", "protein"])
@plac.opt("format_", "return format", choices=["gbwithparts", "fasta", "gb"])
@plac.opt("alias", "remap sequence ids")
def run(db="nuccore", format_="gbwithparts", alias='', *acc):
    ids = []
    for num in acc:
        ids.extend(num.split(","))

    if not sys.stdin.isatty():
        lines = utils.read_lines(sys.stdin, sep=None)
        ids.extend(lines)

    def make_acc(text):
        return Accession(text)

    # Proofread the accession numbers
    accs = map(make_acc, ids)
    for acc in accs:
        if not acc.valid:
            utils.error(f"invalid accession {acc.value}")

        if acc.is_protein():
            db = 'protein'

    ids = ",".join(ids)

    if ids:
        efetch(db=db, rettype=format_, ids=ids)
    else:
        utils.error("no accession numbers were specified")


if __name__ == '__main__':
    # id = "AY851612",

    run()
