import string
import sys, os, io, operator, functools
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from itertools import *
from biorun import utils

NUCS = set("ATGC" + 'atgc')
PEPS = set("ACDEFGHIKLMNPQRSTVWY")
RNAS = set("AUGC" + 'augc')

NUCLEOTIDE, PEPTIDE = "nucleotide", "peptide"

SOURCE = "source"

SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}


def is_nuc(c):
    return c in NUCS


def is_prot(c):
    return c in PEPS


def all_pep(text, limit=1000):
    subs = text[:limit]
    return all(map(is_prot, subs))


def all_nuc(text, limit=1000):
    subs = text[:limit]
    return all(map(is_nuc, subs))


def is_sequence(text):
    if all_nuc(text):
        return NUCLEOTIDE
    elif all_pep(text):
        return PEPTIDE
    else:
        return None


class Peeker:
    def __init__(self, stream):
        self.stream = stream
        self.buffer = io.StringIO()

    def peek(self):
        content = next(self.stream)
        self.buffer = io.StringIO(content)
        return content

    def read(self, size=None):
        if size is None:
            return self.buffer.read() + self.stream.read(size=size)
        content = self.buffer.read(size)
        if len(content) < size:
            content += self.stream.read(size - len(content))
        return content

    def readline(self):
        line = self.buffer.readline()
        if not line.endswith('\n'):
            line += self.stream.readline()
        return line

    def __iter__(self):
        data = self.buffer.read()
        if data:
            yield data
        else:
            for line in self.stream:
                yield line

def get_streams(elems):
    """
    Returns streams including stdin.
    """
    if not sys.stdin.isatty():
        yield sys.stdin

    for fname in elems:
        if os.path.isfile(fname):
            yield open(fname)
        else:
            utils.error(f"file not found: {fname}")

def get_peakable_streams(elems, dynamic=False):
    """
    Returns peekable streams. Can also generate dynamic streams from text.
    """
    label = chain(string.ascii_lowercase)

    if not sys.stdin.isatty():
        yield Peeker(sys.stdin)

    for fname in elems:
        if os.path.isfile(fname):
            yield Peeker(open(fname))
        elif dynamic and is_sequence(fname):
            stream = io.StringIO(f">{next(label)}\n{fname}")
            yield Peeker(stream)
        else:
            utils.error(f"file not found: {fname}")

def parse_stream(stream):
    """
    Guesses the type of input and parses the stream into a BioPython SeqRecord.
    """
    first = stream.peek()
    format = "fasta" if first.startswith(">") else 'genbank'
    recs = SeqIO.parse(stream, format)
    return recs

def flatten(nested):
    return functools.reduce(operator.iconcat, nested, [])


def generate(rec):

    print(type(rec))

    rec.type = SOURCE
    rec.strand = None
    rec.start, rec.end = 1, len(rec.seq)

    yield rec

    #for feat in rec.features:
    #
    #    yield rec




    for feat in rec.features:
        print (type(feat))

def main():

    fnames = sys.argv[1:]

    stream = get_peakable_streams(fnames, dynamic=True)
    reader = map(parse_stream, stream)

    # Flatten all records into a list.
    recs = flatten(reader)

    print (recs)

    recs = map(generate, recs)

    recs = flatten(recs)

    for rec in recs:
        print (rec.type)

if __name__ == '__main__':
    main()
