"""
This package attempts to simplify the BioPython SeqRecords into a simpler, flatter structure that
can be more readily worked with.
"""
import sys, os
from biorun import utils

from itertools import *
from pprint import pprint

try:
    from Bio import SeqIO
except ImportError as exc:
    print(f"*** {exc}", file=sys.stderr)
    print(f"*** This software requires biopython.", file=sys.stderr)
    print(f"*** Try: conda install biopython", file=sys.stderr)
    sys.exit(1)

logger = utils.logger


def get_feature_id(feat):
    """
    Attempts to generate a meaningful name for a feature
    """
    name = utils.first(feat.qualifiers.get("gene", '')) or feat.type

    return name


def get_feature_decription(feat):
    """
    Attempts to generate a meaningful name for a feature
    """
    synon = utils.first(feat.qualifiers.get("gene_synonym", ''))

    db_xref = utils.first(feat.qualifiers.get("db_xref", ''))

    desc = f"type={feat.type} loc={feat.location}"

    if synon:
        desc = f"{desc} synonym={synon}"

    if db_xref:
        desc = f"{desc} db_xref={db_xref}"

    return desc


def filter_features(stream, start=0, end=None, gene=None, typ=None):

    # Filter by type.
    if typ and typ != "all":
        stream = filter(lambda f: f.type == typ, stream)

    # Filter by name.
    if gene:
        stream = filter(lambda f: gene in f.qualifiers.get("gene", []), stream)

    # Filter by start.
    if start:
        stream = filter( lambda f: start <= int(f.location.start), stream)

    # Filter by end.
    if end:
        stream = filter(lambda f: int(f.location.end) <= end, stream)

    return stream


def get_feature_fasta(recs, name='', gene='', start=0, end=None, typ=None):
    """
    Returns records from a list of GenBank
    """

    for item in recs:

        # For a fasta file start and end mean slicing the sequence.
        feats = filter_features(stream=item.features, start=0, end=None, gene=gene, typ=typ)

        for feat in feats:
            # print (feat)
            seq = feat.extract(item)
            seq.id = get_feature_id(feat)
            seq.description = get_feature_decription(feat)
            if start or end:
                end = end or len(item)
                seq.id = f"{seq.id} [{start}:{end}]"
                seq = seq[start: end]
            yield seq


def get_source_fasta(recs, name='', gene='', start=1, end=None, typ=None):
    """
    Returns records from a list of GenBank
    """

    for item in recs:

        # Renane the sequence.
        if name:
            item.id = name

        # If the start/ends are set.
        if start or end:
            end = end or len(item)
            item.id = f"{item.id} [{start}:{end}]"
            item.seq = item.seq[start:end]

        yield item


def parse_genbank(stream, fmt=utils.GENBANK):
    """
    Parses GenBank into sequence records.
    """
    recs = SeqIO.parse(stream, format=fmt)
    return recs


if __name__ == "__main__":
    import doctest

    doctest.testmod()
