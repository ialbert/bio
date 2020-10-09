"""
Fetches data from Entrez.
"""
import sys, os, itertools

import plac
from Bio import Entrez
from biorun import utils
from biorun import models
from biorun.data import fetch
from biorun.utils import first
import itertools

# The default logging function.
logger = utils.logger

# This email needs to be tunable.
Entrez.email = 'bio@bio.com'


def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "mat_peptide": "mature_protein_region",
}

# GFF attributes generated for a source type.
SOURCE_ATTRIBUTES = [
    "mol_type", "isolate", "db_xref", "organism", "country", "collection_date"
]

# GFF attributes filled for each feature other than "source"
GFF_ATTRIBUTES = [
    "gene", "protein_id", "function", "product", "note"
]


def print_fasta(stream, gene='', start=0, end=None, typ=None):
    """
    Prints the origin of a BioPython SeqRecord.
    """
    recs = models.parse_genbank(stream)

    # Print the origin for each record
    if not gene:
        for item in recs:
            print(item.format("fasta"))


def get_pairs(keys, adict):
    pairs = [(k, adict.get(k)) for k in keys]
    pairs = [f"{k}={first(v)}" for (k, v) in pairs if v]
    return pairs


def make_attr(feat):
    """
    Creates GFF attributes from a SeqRecord
    """
    gene = first(feat.qualifiers.get("gene"))
    data = []

    if feat.type == "gene":
        data.append(f"Name={gene}")

    if feat.type == "source":
        data.extend(get_pairs(SOURCE_ATTRIBUTES, feat.qualifiers))
    else:
        # Adds generic attrbiutes
        data.extend(get_pairs(GFF_ATTRIBUTES, feat.qualifiers))

    return ";".join(data)


def feature2gff(feat, anchor):
    """
    Returns a SeqRecord as an 11 element  GFF3 list .
    """
    start = int(feat.location.start)
    end = int(feat.location.end)
    attr = make_attr(feat)
    strand = "+" if feat.strand else "-"
    ftype = SEQUENCE_ONTOLOGY.get(feat.type, feat.type)
    data = [anchor, ".", ftype, start, end, ".", strand, ".", attr]
    return data


def filter_features(stream, start=0, end=None, gene=None, typ=None):
    # Filter by type.
    if typ:
        def type_filter(f):
            return f.type == typ

        stream = filter(type_filter, stream)

    # Filter by name.
    if gene:
        def name_filter(f):
            return gene in f.qualifiers.get("gene", [])

        stream = filter(name_filter, stream)

    # Filter by start.
    if start:
        def start_filter(f):
            return start <= int(f.location.start)

        stream = filter(start_filter, stream)

    # Filter by end.
    if end:
        def end_filter(f):
            return int(f.location.end) <= end

        stream = filter(end_filter, stream)

    return stream


def print_gff(stream, gene='', start=0, end=None, typ=None):
    """
    Prints the origin of a BioPython SeqRecord.
    """
    recs = models.parse_genbank(stream)

    # Print the origin for each record
    for rec in recs:

        # Subselect by coordinates.
        feats = filter_features(stream=rec.features, start=start, end=end, gene=gene, typ=typ)

        # Generate the gff output
        for feat in feats:
            values = feature2gff(feat, anchor=rec.id)
            values = map(str, values)
            print("\t".join(values))


def print_genbank(stream):
    for line in stream:
        print(line, end="")


def process(acc, gene='', fasta=False, gff=False, start=0, end=None, typ=''):
    """
    Performs the processing of a single accession number.
    """

    # Open the stream to the data
    cache, stream = fetch.get(acc=acc)

    if fasta:
        print_fasta(stream, gene=gene, start=start, end=end, typ=type)
    elif gff:
        print_gff(stream, gene=gene, start=start, end=end, typ=typ)
    else:
        print_genbank(stream)

    return


@plac.pos('acc', "accession numbers")
@plac.opt('gene', "name of the gene associated with the feature")
@plac.opt('type', "the type of the feature")
@plac.opt('start', "start coordinate ")
@plac.opt('end', "end coordinate")
@plac.flg('fasta', "generate fasta file")
@plac.flg('gff', "generate a gff file", abbrev="G")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(gene='', type='', start=0, end=0, fasta=False, gff=False, verbose=False, *acc):
    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    start = start
    end = end or None
    # Process each accession number.
    for acx in acc:
        process(acx, gene=gene, fasta=fasta, gff=gff, start=start, end=end, typ=type)
