"""
This package attempts to simplify the BioPython SeqRecords into a simpler, flatter structure that
can be more readily worked with.
"""
import sys, os, json
from biorun import utils
from collections import OrderedDict

from itertools import *
from pprint import pprint

from biorun.const import *



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


def filter_features(items, start=0, end=None, gene=None, typ=None, translation=None):
    # Filter by type.
    if typ and typ != "all":
        items = filter(lambda f: f.get('type') == typ, items)

    # Filter by name.
    if gene:
        items = filter(lambda f: gene in f.get("gene", []), items)

    # Filter by the translation attribute existance.
    if translation:
        items = filter(lambda f: f.get("translation", [''])[0], items)

    # Filter by start.
    if start:
        items = filter(lambda f: start <= f.get('start'), items)

    # Filter by end.
    if end:
        items = filter(lambda f: f.get('end') <= end, items)

    return items


def get_feature_fasta(data,  gene='', start=0, end=None, typ=None, translation=None):
    """
    Returns records from a list of GenBank
    """

    elems = data[FEATURES]
    elems = filter_features(elems, start=0, end=None, gene=gene, typ=typ, translation=translation)

    for elem in elems:
        if translation:
            text = elem.get("translation", [''])[start:end]
            rec = SeqRecord(Seq(text), id="foo")
            yield rec


    '''
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
    '''


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

def get_origin(item, start=0, end=None, name=None):
    """
    Returns the origin sequence from an JSON item
    """
    # Prints the source sequence
    seq = item[ORIGIN][start:end]
    desc = item[DEFINITION]
    seqid = item[SEQID]
    locus = item[LOCUS]
    seqid = name or seqid
    rec = SeqRecord(Seq(seq), id=seqid, name=locus, description=desc)
    return rec

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
        item[SEQID] = rec.id
        item[DEFINITION] = rec.description
        item[DBLINK] = rec.dbxrefs
        item[LOCUS] = rec.name

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
        item[ORIGIN] = str(rec.seq)

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
