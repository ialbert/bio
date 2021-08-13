import sys, re

try:
    from Bio import Entrez
except ImportError as exc:
    print(f"### Error: {exc}", file=sys.stderr)
    print(f"### This program requires biopython", file=sys.stderr)
    print(f"### Install: conda install -y biopython>=1.79", file=sys.stderr)
    sys.exit(-1)

from biorun.libs import placlib as plac
from tqdm import tqdm
from biorun import utils

Entrez.email = 'foo@foo.com'

ncbi_patt = r'(?P<letters>[a-zA-Z]+)(?P<under>_?)(?P<digits>\d+)(\.(?P<version>\d+))?'
ncbi_patt = re.compile(ncbi_patt)

def detect_format(text):

    # Allow integer (gi numbers)
    try:
        int(text)
        return text
    except ValueError as exc:
        pass

    m = ncbi_patt.search(text)
    if not m:
        utils.error(f"accession number format not recognized: {text}")

    # Unused at this time.
    letters, digits, under, version = m.group("letters"), m.group("digits"), m.group("under"), m.group("version")

    return text


def efetch(ids, db='nuccore', rettype='gbwithparts', retmode='text'):
    stream = Entrez.efetch(db="nucleotide", id=ids, rettype=rettype, retmode=retmode)
    stream = tqdm(stream, unit='B', unit_divisor=1024, desc='# downloaded', unit_scale=True, delay=5, leave=False)
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
        lines = utils.read_lines(sys.stdin, sep='')
        ids.extend(lines)

    ids = map(detect_format, ids)
    ids = ",".join(ids)

    if ids:
        efetch(db=db, rettype=format_, ids=ids)
    else:
        utils.error("no accession numbers were specified")

if __name__ == '__main__':
    # id = "AY851612",

    run()
