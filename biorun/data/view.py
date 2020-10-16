"""
Fetches data from Entrez.
"""
import sys, os, itertools, json, re

import plac

from biorun import models
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from biorun import utils
from biorun.data import fetch
from biorun.utils import first
from biorun.const import *

import itertools

# The default logging function.
logger = utils.logger

# This email needs to be tunable.
Entrez.email = 'bio@bio.com'


def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


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

def get_pairs(keys, adict):
    pairs = [(k, adict.get_data(k)) for k in keys]
    pairs = [f"{k}={first(v)}" for (k, v) in pairs if v]
    return pairs


def make_attr(feat):
    """
    Creates GFF attributes from a SeqRecord
    """
    gene = first(feat.qualifiers.get_data("gene"))
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
    ftype = SEQUENCE_ONTOLOGY.first(feat.type, feat.type)
    data = [anchor, ".", ftype, start, end, ".", strand, ".", attr]
    return data


def print_gff(stream, gene='', name='', start=0, end=None, typ=None):
    """
    Prints the origin of a BioPython SeqRecord.
    """
    recs = models.parse_genbank(stream)

    # Print the origin for each record
    for rec in recs:

        # The name of the GFF anchor.
        anchor = name or rec.id

        # Subselect by coordinates.
        feats = models.filter_features(stream=rec.features, start=start, end=end, gene=gene, ftype=typ)

        # Generate the gff output
        for feat in feats:
            values = feature2gff(feat, anchor=anchor)
            values = map(str, values)
            print("\t".join(values))


def process(acc, param):
    """
    Performs the processing of a single accession number.
    """

    # Open the stream to the data
    data = fetch.get_data(acc=acc)

    # When to produce the origin fasta.
    origin = param.fasta and not(param.gene or param.type or param.protein or param.translate)

    if param.protein:
        print_translation_fasta(data, param=param)
    elif origin:
        print_origin_fasta(data, param=param)
    elif param.fasta:
        print_feature_fasta(data, param=param)
    elif param.gff:
        print_gff(data, gene=gene, seqid=seqid, start=start, end=end, typ=typ)
    else:
        print_json(data, param=param)


    return


@plac.pos('acc', "accession numbers")
@plac.opt('gene', "name of the gene associated with the feature")
@plac.opt('seqid', "set the sequence id", abbrev="Q")
@plac.opt('match', "select elements by matching a regexp to any existing information")
@plac.opt('type', "filter by the type of the feature")
@plac.opt('start', "start coordinate ", type=int)
@plac.opt('end', "end coordinate", type=int)
@plac.flg('fasta', "generate fasta file")
@plac.flg('protein', "fasta file with protein translations embedded in the data")
@plac.flg('translate', "translates DNA sequences to protein", abbrev="T")
@plac.flg('gff', "generate a gff file", abbrev="G")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(gene='', type='', seqid='', match='', start=1, end=0, fasta=False, protein=False, translate=False, gff=False, verbose=False, *acc):
    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    if start == 0:
        error(f"start={start} may not be  zero")

    # Move it to one based coordinate system
    start = start - 1 if start > 0 else start

    # Set the end slice
    end = end or None

    # Set the regular expression.
    regexp = re.compile(match) if match else match

    # Collected parameters
    param = utils.Param(start=start, end=end, seqid=seqid, protein=protein,
                  gff=gff, translate=translate, fasta=fasta, type=type, gene=gene, regexp=regexp)

    # Process each accession number.
    for acx in acc:
        process(acx, param=param)
