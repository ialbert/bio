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

def print_origin_fasta(data,  param):
    """
    Prints the origin sequence for each record.
    """
    for item in data:
        rec = models.get_origin(item, param)
        print(rec.format("fasta"))

def print_translation_fasta(data, param):
    """
    Prints the translation field of each record.
    """
    for item in data:
        recs = models.get_translation_records(item, param)
        for rec in recs:
            print(rec.format("fasta"))

def print_feature_fasta(data,  param):
    """
    Prints the sequence for features.
    """

    for item in data:
        recs = models.get_feature_records(item, param)
        for rec in recs:
            print(rec.format("fasta"))

def print_json(data,  param):
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
        feats = models.filter_features(feats, start=param.start, end=param.end, ftype=param.type, gene=param.gene, regexp=param.regexp)
        text = json.dumps(list(feats), indent=4)
        print (text)


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

    print ("##gff-version 3")
    
    for item in data:

        feats = item[FEATURES]

        # The name of the GFF anchor.
        anchor = param.seqid or item['id']

        # Subselect by coordinates.
        feats = models.filter_features(feats, start=param.start, end=param.end, gene=param.gene, ftype=param.type, regexp=param.regexp)

        # Generate the gff output
        for feat in feats:
            values = feature2gff(feat, anchor=anchor)
            values = map(str, values)
            print("\t".join(values))


def convert_all(names, param):
    """
    Converts all accession numbers listed
    """
    for name in names:
        data = storage.get_json(name, seqid=param.seqid)
        if not data:
            utils.error(f"data not found: {name}")
        convert_one(data, param)

def convert_one(data, param):
    """
    Converts an accession number
    """

    # When to produce the origin fasta.
    origin = param.fasta and not(param.gene or param.type or param.protein or param.translate)

    if param.protein:
        print_translation_fasta(data, param=param)
    elif origin:
        print_origin_fasta(data, param=param)
    elif param.gff:
        print_gff(data, param=param)
    elif param.fasta or param.translate:
        print_feature_fasta(data, param=param)
    else:
        print_json(data, param=param)

