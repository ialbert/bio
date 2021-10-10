import functools
import io
import operator
import os
import sys
import gzip
import json
from pprint import pprint

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import Reference, CompoundLocation, FeatureLocation
except ImportError as exc:
    print(f"# Error: {exc}", file=sys.stderr)
    print(f"# This program requires biopython", file=sys.stderr)
    print(f"# Install: conda install -y biopython>=1.79", file=sys.stderr)
    sys.exit(-1)

from collections import OrderedDict, defaultdict
from itertools import *
from biorun import utils

logger = utils.logger

NUCS = set("ATGC" + 'atgc')
PEPS = set("ACDEFGHIKLMNPQRSTVWY")
RNAS = set("AUGC" + 'augc')

NUCLEOTIDE, PEPTIDE = "nucleotide", "peptide"

SOURCE, FEATURES, RECORD, ID = "source", "features", "record", "id",
TYPE, ANNOTATIONS, LOCATIONS =  "type", "annotations", "locations"
SEQUENCE = "sequence"

SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}

counter = count(1)

UNIQUE = defaultdict(int)


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
        try:
            content = next(self.stream)
        except StopIteration:
            content = ''
        self.buffer = io.StringIO(content)
        return content

    def read(self, size=None):
        if size is None:
            return self.buffer.read() + self.stream.read()
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

    def __next__(self):
        data = self.buffer.read() or next(self.stream)
        if not data:
            raise StopIteration
        else:
            return data


def get_streams(fnames, dynamic=False):
    """
    Yields streams. stdin is generated first
    """
    label = count(1)

    if not sys.stdin.isatty():
        logger.debug(f"reading stdin")
        yield sys.stdin

    for fname in fnames:
        fname = fname.strip()
        if os.path.isfile(fname):
            if fname.endswith(".gz"):
                logger.debug(f"gzip open: {fname}")
                yield gzip.open(fname, mode='rt')
            else:
                logger.debug(f"open: {fname}")
                yield open(fname)
        elif dynamic and is_sequence(fname):
            logger.debug(f"dynamic sequence")
            stream = io.StringIO(f">seq{next(label)}\n{fname}")
            yield stream
        else:
            if fname.startswith("--"):
                utils.error(f"invalid parameter: {fname}")
            else:
                utils.error(f"file not found: '{fname}'")

def parse_json(stream):
    text = stream.read()
    data = json.loads(text)
    for entry in data:
        feats = entry[FEATURES]
        parents = [f for f in feats if f[TYPE] == SOURCE]
        if parents:
            parent_ann = parents[0][ANNOTATIONS]
        else:
            parent_ann = {}

        for feat in entry[FEATURES]:
            uid = feat[ID]
            ftype = feat[TYPE]
            ann = feat[ANNOTATIONS]
            desc = ""
            seq = Seq("ATGC")
            rec = BioRec(id=uid, ann=ann, parent=parent_ann, type=ftype, seq=seq, desc=desc)
            yield rec

def parse_stream(stream):
    """
    Guesses the type of input and parses the stream into a BioPython SeqRecord.
    """
    stream = Peeker(stream)
    first = stream.peek().strip()
    start = first[0] if first else ''
    if start == '>':
        format = 'fasta'
    elif start == '@':
        format = 'fastq'
    elif first.startswith("ID "):
        format = 'embl'
    elif start == '[':
        format = 'json'
    else:
        format = 'genbank'

    logger.debug(f"parsing: {format}")

    if format != 'json':
        recs = SeqIO.parse(stream, format)
    else:
        recs = parse_json(stream)

    return recs


def flatten(nested):
    return functools.reduce(operator.iconcat, nested, [])


def json_ready(value):
    """
    Recursively serializes values to types that can be turned into JSON.
    """

    # The type of of the incoming value.
    curr_type = type(value)

    # Reference types will be dictionaries.
    if curr_type == Reference:
        return dict(title=value.title, authors=value.authors, journal=value.journal, pubmed_id=value.pubmed_id)

    # Serialize each element of the list.
    if curr_type == list:
        return [json_ready(x) for x in value]

    # Serialize the values of an Ordered dictionary.
    if curr_type == OrderedDict:
        return dict((k, json_ready(v)) for (k, v) in value.items())

    return value


def next_count(ftype):
    UNIQUE[ftype] += 1
    return f'{ftype}-{UNIQUE[ftype]}'


def first(data, key, default=""):
    # First element of a list value that is stored in a dictionary by a key.
    value = data.get(key)
    if type(value) == list:
        return value[0]
    elif value:
        return value
    else:
        return default


def guess_name(ftype, annot):
    """
    Attempts to generate an unique id name for a BioPython feature
    """
    uid = ''

    gene = first(annot, "gene")
    prod = first(annot, "product")
    locus = first(annot, "locus_tag")
    data = dict(type=ftype, gene=gene, product=prod, locus=locus)

    if ftype == 'gene':
        name = first(annot, "gene") or locus

    elif ftype == 'CDS':
        name = first(annot, "protein_id")
    elif ftype == 'mRNA':
        name = first(annot, "transcript_id")
    elif ftype == "exon":
        name = gene = first(annot, "gene")
        uid = f"{gene}-{ftype}"
        data['number'] = first(annot, "number") or 1
    else:
        name = next_count(ftype)

    desc = json.dumps(data)
    name = name or next_count(f"id")
    uid = uid or name
    return uid, name, desc


def get_id(ftype, annot):
    """
    Attempts to generate an unique id name for a BioPython feature
    """

    gene = first(annot, "gene")
    prod = first(annot, "product")
    locus = first(annot, "locus_tag")

    data = dict(type=ftype, gene=gene, product=prod, locus=locus)

    if ftype == 'gene':
        name = first(annot, "gene") or locus
    elif ftype == 'CDS':
        name = first(annot, "protein_id")
    elif ftype == 'mRNA':
        name = first(annot, "transcript_id")
    else:
        name = next_count(ftype)

    desc = json.dumps(data)

    return name, desc


class BioRec:

    def __init__(self, id, type=None, source=None, ann=None, parent=None, seq=None, desc=''):
        self.id = self.name = id
        self.type = type or SOURCE
        self.source = source or self.id
        self.strand = 0
        self.parent = parent or {}
        self.ann = ann or {}
        self.seq = seq or Seq('')
        self.start, self.end = 1, len(self.seq)
        self.locs = []
        self.gene = first(self.ann, "gene")
        self.product = first(self.ann, "product")
        self.desc = desc or ''

    def get_ann(self, name):
        return first(self.ann, name)

    def get_parent_ann(self, name):
        return first(self.parent, name)

    def __repr__(self):
        return f"{self.__class__.__name__}(id={self.id})"

    def __len__(self):
        return len(self.seq)

def parse_desc(srec):
    """
    Parses annotations from a sequence description
    """
    ann = {}
    desc = srec.description

    # Attempts to split the description for parsing into JSON
    payload = srec.description.split(" ", maxsplit=1)
    if len(payload) == 2:
        try:
            desc = payload[1]
            ann = json.loads(desc)
        except Exception:
            # Not a valid
            pass

    return ann, desc


def record_generator(rec):
    """
    Creates BioRec from SeqRec
    """

    # Already in the correct format
    if isinstance(rec, BioRec):
        yield rec
        raise StopIteration

    # Handling single or nested records.
    if not rec.features:
        # Attempts to read annotations from description.
        ann, desc = parse_desc(rec)
        brec = BioRec(id=rec.id, ann=ann, parent=ann, type=SOURCE, seq=rec.seq, desc=desc)
        yield brec

    else:
        # Reads annotations from the sequence record
        pairs = [(k, json_ready(v)) for (k, v) in rec.annotations.items()]

        # These are the parent annotations that apply globally to all features.
        parent_ann = dict(pairs)

        # Parent record
        parent_rec = rec

        # Yield a record for each feature.
        for feat in rec.features:

            # Extract feature qualifiers.
            pairs = [(k, json_ready(v)) for (k, v) in feat.qualifiers.items()]
            ann = dict(pairs)

            # Remap types into standard naming.
            feat.type = SEQUENCE_ONTOLOGY.get(feat.type, feat.type)

            # Extract the sequence for the feature.
            seq = feat.extract(rec.seq)

            # The source annotations are global
            if feat.type == SOURCE:
                parent_ann.update(ann)
                uid = parent_rec.id
                #isolate = first(parent_ann, "isolate")
                data = dict(title=parent_rec.description, type=SOURCE)
                desc = json.dumps(data)
            else:
                # Find the unique id of the feature.
                uid, desc = get_id(feat.type, ann)

            # Create the BioRecord for the feature.
            brec = BioRec(id=uid, ann=ann, parent=parent_ann, seq=seq, type=feat.type, desc=desc, source=parent_rec.id)

            # Correct the feature coordinates.
            brec.strand = feat.strand

            brec.start, brec.end = int(feat.location.start) + 1, int(feat.location.end)
            brec.locs = [(loc.start + 1, loc.end, loc.strand) for loc in
                        feat.location.parts]

            yield brec


def rec2seqrec(rec):
    """
    Creates Seqrecord from BioRec
    """
    srec = SeqRecord(id=rec.id, name=rec.id, description=rec.desc, seq=rec.seq)
    return srec


def get_records(fnames):
    """
    Create a single stream of SeqRecords from multiple sources
    """
    stream = get_streams(fnames, dynamic=True)

    reader = map(parse_stream, stream)

    recs = flatten(reader)

    recs = map(record_generator, recs)

    recs = flatten(recs)

    return recs


def main():
    fnames = sys.argv[1:]

    stream = get_streams(fnames, dynamic=True)
    reader = map(parse_stream, stream)
    recs = flatten(reader)

    recs = map(record_generator, recs)
    recs = flatten(recs)

    for rec in recs:
        print(rec.name, rec.type)


if __name__ == '__main__':
    main()
