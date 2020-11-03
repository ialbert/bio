"""
Fetches data from Entrez.
"""
import sys, json

from biorun import models
from biorun import utils
from biorun.data import storage
from biorun.const import *

# The default logging function.
logger = utils.logger


def get_origin_fasta(data, param):
    """
    Prints the origin sequence for each record.
    """
    for item in data:
        rec = models.get_origin(item, param)
        yield rec


def get_translation_fasta(data, param):
    """
    Prints the translation field of each record.
    """
    for item in data:
        recs = models.get_translation_records(item, param)
        for rec in recs:
            yield rec


def get_feature_fasta(data, param):
    """
    Prints the sequence for features.
    """
    for item in data:
        recs = models.get_feature_records(item, param)
        for rec in recs:
            yield rec

def get_fasta(data, param):
    """
    Returns SeqRecord objects based on parameters.
    """
    # If there is no other filtering, produce the origin.
    origin = not (param.gene or param.type or param.protein or param.translate)

    if origin:
        recs = get_origin_fasta(data, param=param)
    elif param.protein:
        recs = get_translation_fasta(data, param=param)
    else:
        recs = get_feature_fasta(data, param=param)

    return recs

def print_fasta(recs):
    """
    Prints fasta records
    """
    for rec in recs:
        print(rec.format("fasta"))


def print_json(data, param):
    """
    Prints the sequence for features.
    """

    # Produce the full file when no parameters are set.
    if param.unset():
        text = json.dumps(data, indent=4)
        print(text)
        return

    # Selects individual features.
    for item in data:
        feats = item[FEATURES]
        feats = models.filter_features(feats, start=param.start, end=param.end, ftype=param.type, gene=param.gene,
                                       regexp=param.regexp)
        text = json.dumps(list(feats), indent=4)
        print(text)


def feature2gff(feat, anchor):
    """
    Returns a SeqRecord as an 11 element  GFF3 list .
    """
    start = feat['start']
    end = feat['end']
    ftype = feat['type']
    strand = feat['strand']
    phase = feat.get("codon_start", [1])[0]
    attr = models.make_attr(feat)
    strand = "+" if strand else "-"
    ftype = SEQUENCE_ONTOLOGY.get(ftype, ftype)
    data = [anchor, ".", ftype, start, end, ".", strand, phase, attr]
    return data


def print_gff(data, param):
    """
    Prints the origin of a BioPython SeqRecord.
    """

    print("##gff-version 3")

    for item in data:

        feats = item[FEATURES]

        # The name of the GFF anchor.
        anchor = param.seqid or item['id']

        # Subselect by coordinates.
        feats = models.filter_features(feats, start=param.start, end=param.end, gene=param.gene, ftype=param.type,
                                       regexp=param.regexp)

        # Generate the gff output
        for feat in feats:
            values = feature2gff(feat, anchor=anchor)
            values = map(str, values)
            print("\t".join(values))


def convert(names, params):
    """
    Converts all accession numbers listed
    """
    for name, param in zip(names, params):
        data = storage.get_json(name, seqid=param.seqid)
        if not data:
            utils.error(f"data not found: {name}")
        convert_one(data, param)


def convert_one(data, param):
    """
    Converts an accession number
    """

    # GFF conversion.
    if param.gff:
        print_gff(data, param=param)
        return

    # Some type of FASTA conversion
    fasta = param.fasta or param.protein or param.translate
    if fasta:
        recs = get_fasta(data, param=param)
        print_fasta(recs)
        return

    # No explicit conversion request, print JSON
    print_json(data, param=param)
