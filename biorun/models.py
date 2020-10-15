"""
This package attempts to simplify the BioPython SeqRecords into a simpler, flatter structure that
can be more readily worked with.
"""
import sys, os, json
from biorun import utils
from collections import OrderedDict

from itertools import *
from pprint import pprint

# The key by under which the features are stored
FEATURES = "FEATURES"

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import Reference, CompoundLocation, FeatureLocation
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
    name = utils.first(feat.qualifiers.get_data("gene", '')) or feat.type

    return name


def get_feature_decription(feat):
    """
    Attempts to generate a meaningful name for a feature
    """
    synon = utils.first(feat.qualifiers.get_data("gene_synonym", ''))

    db_xref = utils.first(feat.qualifiers.get_data("db_xref", ''))

    desc = f"type={feat.type} loc={feat.location}"

    if synon:
        desc = f"{desc} synonym={synon}"

    if db_xref:
        desc = f"{desc} db_xref={db_xref}"

    return desc


def filter_features(stream, start=0, end=None, gene=None, typ=None, translation=None):
    # Filter by type.
    if typ and typ != "all":
        stream = filter(lambda f: f.type == typ, stream)

    # Filter by name.
    if gene:
        stream = filter(lambda f: gene in f.qualifiers.get_data("gene", []), stream)

    # Filter by the translation attribute existance.
    if translation:
        stream = filter(lambda f: f.qualifiers.get_data("translation"), stream)

    # Filter by start.
    if start:
        stream = filter(lambda f: start <= int(f.location.start), stream)

    # Filter by end.
    if end:
        stream = filter(lambda f: int(f.location.end) <= end, stream)

    return stream


def get_feature_fasta(recs, name='', gene='', start=0, end=None, typ=None, translation=None):
    """
    Returns records from a list of GenBank
    """

    for item in recs:

        # For a fasta file start and end mean slicing the sequence.
        feats = filter_features(stream=item.features, start=0, end=None, gene=gene, typ=typ, translation=translation)

        for feat in feats:
            # print (feat)
            if translation:
                text = utils.first(feat.qualifiers.get_data("translation", ''))
                seq = SeqRecord(Seq(text), id="seq")
            else:
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


def serialize(value):
    """
    Serializes values to JSON ready type.
    """

    # The type of of the incoming value.
    curr_type = type(value)

    # Reference types will be dictionaries.
    if curr_type == Reference:
        return dict(title=value.title, authors=value.authors, journal=value.journal, pubmed_id=value.pubmed_id)

    # Serialize each element of the list.
    if curr_type == list:
        return [serialize(x) for x in value]

    # Serialize the elements of an Ordered dictionary.
    if curr_type == OrderedDict:
        return dict((k, serialize(v)) for (k, v) in value.items())

    return value


@utils.time_it
def convert_genbank(stream):
    """
    Converts a stream to a GenBank file into json.
    """

    # Parses the genbank into memory.
    recs = parse_genbank(stream)

    # The outer dictionary containing multiple records.
    data = []

    # Add each record separately.
    for rec in recs:

        # The individual SeqRecords as dictionaries.
        item = dict()

        # Fill the standard SeqRecord fields.
        item['id'] = rec.id
        item['definition'] = rec.description
        item['dblink'] = rec.dbxrefs
        item['locus'] = rec.name

        # Fill in all annotations.
        for key, value in rec.annotations.items():
            item[key] = serialize(value)

        # Fill in the features.
        feats = []
        for feat in rec.features:
            ftype = feat.type
            start = int(feat.location.start) + 1
            end = int(feat.location.end)

            # print (feat.location, type(feat.location))
            location = [(loc.start + 1, loc.end, loc.strand) for loc in feat.location.parts]
            elem = dict(start=start, end=end, type=ftype, location=location)

            for (k, v) in feat.qualifiers.items():
                elem[k] = serialize(v)

            feats.append(elem)

        # Add the features
        item[FEATURES] = feats

        # Save the sequence as well
        item['seq'] = str(rec.seq)

        # Each item keyed as the record.id.
        data.append(item)

        # pprint(rec.annotations)

        # print(rec)

    return data


def parse_genbank(stream, fmt=utils.GENBANK):
    """
    Parses GenBank into sequence records.
    """
    recs = SeqIO.parse(stream, format=fmt)
    return recs


if __name__ == "__main__":
    import doctest

    doctest.testmod()
