"""
Fetches data from Entrez.
"""
import sys, os, itertools, json

import plac

from biorun import models
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from biorun import utils
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


def print_fasta(data, name='', gene='', start=0, end=None, typ=None, translation=None):
    """
    Prints the origin of a BioPython SeqRecord.
    """

    for item in data:
        if gene or typ or translation:
            seqs = models.get_feature_fasta(item, start=start, end=end, name=name, gene=gene, typ=typ, translation=translation)
        else:
            seq = item['seq'][start:end]
            seqid = item['id']
            desc = item['definition']
            name = name or seqid
            seq = SeqRecord(Seq(seq), id=seqid, name=name, description=desc)
            seqs = [ seq ]

        # Print the sequence for each record
        for elem in seqs:
            print(elem.format("fasta"))


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
    ftype = SEQUENCE_ONTOLOGY.get(feat.type, feat.type)
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
        feats = models.filter_features(stream=rec.features, start=start, end=end, gene=gene, typ=typ)

        # Generate the gff output
        for feat in feats:
            values = feature2gff(feat, anchor=anchor)
            values = map(str, values)
            print("\t".join(values))


def print_genbank(stream):
    for line in stream:
        print(line, end="")


def process(acc, gene='', name='', fasta=False, gff=False, start=0, end=None, typ='', translation=None):
    """
    Performs the processing of a single accession number.
    """

    # Open the stream to the data
    data = fetch.get_data(acc=acc)

    # Turns on fasta formats
    fasta = fasta or (not gff and gene)

    if fasta or translation:
        print_fasta(data, gene=gene,  name=name, start=start, end=end, typ=typ, translation=translation)
    elif gff:
        print_gff(data, gene=gene, name=name, start=start, end=end, typ=typ)
    else:
        text = json.dumps(data, indent=4)
        print(text)

    return


@plac.pos('acc', "accession numbers")
@plac.opt('gene', "name of the gene associated with the feature")
@plac.opt('rename', "renames the elements")
@plac.opt('type', "the type of the feature")
@plac.opt('start', "start coordinate ", type=int)
@plac.opt('end', "end coordinate", type=int)
@plac.flg('fasta', "generate fasta file")
@plac.flg('translation', "get the features that have a translation attribute set", abbrev="T")
@plac.flg('gff', "generate a gff file", abbrev="G")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(gene='', type='', rename='', start=1, end=None, fasta=False, translation=False, gff=False, verbose=False, *acc):
    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    if start < 1:
        error(f"start={start} must be greater than zero")

    if (end is not None) and (start > end):
        error(f"start={start} is larger than end={end}")

    # Move it to one based coordinate system
    start = start - 1

    # Set thethe end slice
    end = end or None

    # Process each accession number.
    for acx in acc:
        process(acx, gene=gene, name=rename, fasta=fasta, gff=gff, start=start, end=end, typ=type, translation=translation)
