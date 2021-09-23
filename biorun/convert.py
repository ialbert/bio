"""
Converts across formats
"""
import gzip
import os
import sys
from collections import OrderedDict, defaultdict
from itertools import count, tee
from Bio.SeqIO.FastaIO import FastaIterator

from biorun import utils, parser

# Module level logger.
logger = utils.logger

# Global alias remapper
ALIAS = {}

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import Reference, CompoundLocation, FeatureLocation
except ImportError as exc:
    logger.error(f"{__name__}", stop=False)
    logger.error(f"{exc}", stop=False)
    logger.error(f"This software requires biopython.")
    logger.error(f"Try: conda install biopython>=1.79")
    sys.exit()


def fasta_formatter(rec):
    print(rec.format("fasta"), end='')


def remapper(rec):
    rec.id = ALIAS.get(rec.id, rec.id)
    return rec


def interval_selector(start=0, end=None):
    def func(rec):
        if start and rec.start < start:
            return False
        if end and rec.end > end:
            return False
        return True

    return func


def sequence_slicer(start=0, end=None):
    def func(rec):
        if start or end:
            seqlen = len(rec.seq)

            endx = seqlen if end is None else end
            endx = seqlen + endx if endx < 0 else endx

            # Zero based shift
            if start < 0:
                startx = seqlen + start + 1
            else:
                startx = start + 1

            rec.description = f"{rec.description} [{startx}:{endx}]"
            rec.seq = rec.seq[start:end]
        return rec

    return func

def true_all(x):
    return True

def type_selector(ftype):

    if ftype =='all':
        return lambda x: True

    types = set(ftype.split(","))

    def func(rec):
        return rec.type in types if ftype else True

    return func


def gene_selector(name):
    genes = set(name.split(","))

    def func(rec):
        return parser.first(rec.annot, "gene") in genes if name else True

    return func


def name_selector(name):
    names = set(name.split(","))

    def func(rec):
        return rec.id in names if name else True

    return func


def seqid_selector(seqid):
    targets = set(seqid.split(","))

    def func(rec):
        return rec.id in targets if seqid else True

    return func


def translate_recs(flag):
    def func(rec):
        if flag:
            endx = len(rec.seq) // 3
            rec.seq = rec.seq[0:endx * 3].translate()
        return rec

    return func


def reverse_complement(flag):
    """
    Reverse complement
    """

    def func(rec):
        if flag:
            rec.seq = rec.seq.reverse_complement()
        return rec

    return func


def source_only(flag):
    def func(rec):
        return rec.type == parser.SOURCE if flag else True

    return func

def keep_source(rec):
    return rec.type == parser.SOURCE


def remove_source(rec):
    return rec.type != parser.SOURCE


def protein_filter(rec):
    return "translation" in rec.annot


def protein_extract(rec):
    rec.seq = Seq(rec.annot.get("translation")[0])
    return rec


#
# The GFF attributes generated for a source type.
#
SOURCE_ATTRIBUTES = [
    "mol_type", "isolate", "db_xref", "organism", "country", "collection_date"
]

# GFF attributes filled for each feature other than "source"
GFF_ATTRIBUTES = [
    "gene", "protein_id", "product", "db_xref", "function",
]

# Attributes that should not be added
SKIP_GFF_ATTR = {"id", "parent_id", "name", "type", "start", "comment", "references", "structured_comment",
                 "end", "location", "translation", "strand", "operator"}

# Associate a color to a feature type.
COLOR_FOR_TYPE = {
    "five_prime_UTR": "#cc0e74",
    "three_prime_UTR": "#cc0e74",
    "stem_loop": "#fa7f72",
    "mature_protein_region": "#CBAEBB",
    "region": "#CECECE",
    "mRNA": "#799351",
    "gene": "#cb7a77",
    "transcript": "#79a3b1",
    "tRNA": "#a685e2",
    "ncRNA": "#fca3cc",
    "mobile_element": "#efd9d1",
    "mRNA_region": "#7a77cb",
}


def feature2gff(anchor, ftype, start, end, strand, uid, name, pid=None):
    """
    Returns a Record as an 11 element  GFF3 list .
    """
    # Reformat the strand
    strand = "+" if strand > 0 else "-"

    # TODO: is this the phase?
    # phase = feat.get("codon_start", [1])[0] - 1
    phase = "."

    # The color for the feature.
    color = COLOR_FOR_TYPE.get(ftype)

    # Attribute data
    attr = [f"ID={uid}", f"Name={name}"]
    if pid:
        attr.append(f"Parent={pid}")
    if color:
        attr.append(f"color={color}")

    # Build the attribute string
    attr = ";".join(attr)

    # Create the GFF record.
    data = [anchor, ".", ftype, start, end, ".", strand, phase, attr]

    return data


def gff_formatter(rec):
    """
    Formats a record as GFF.
    """

    if rec.type == parser.SOURCE:
        return

    # Parent feature
    data = feature2gff(start=rec.start, end=rec.end,
                       ftype=rec.type, uid=rec.id,
                       name=rec.name, strand=rec.strand,
                       anchor=rec.anchor, pid=None)

    line = "\t".join(map(str, data))

    # Parent id.
    pid = rec.id

    if rec.type == "mRNA":
        print(line)
        ftype = "exon"
    else:
        ftype = rec.type

    # Generate the locations
    for start, end, strand in rec.locs:
        name = rec.name
        uid = next(parser.counter)

        data = feature2gff(start=start, end=end, ftype=ftype, uid=uid, name=name,
                           strand=strand, anchor=rec.anchor, pid=rec.id)
        line = "\t".join(map(str, data))
        print(line)


def run(features=False, protein=False, translate=False, fasta=False, revcomp=False,
        start='1', end=None, type_='', id_='', name='', gene='', alias=None, fnames=[]):
    """
    Converts data to different formats.
    """
    global ALIAS

    # Generate the ALIAS remapping.
    ALIAS = utils.parse_alias(alias) if alias else {}

    # Parse start and end into user friendly numbers.
    start = utils.parse_number(start)

    # Positive coordinates moved to 0 based.
    if start >= 0:
        start = start - 1

    # End coordinate.
    end = utils.parse_number(end)

    # Set additional filtering paramters.
    ftype = type_
    seqid = id_

    # Turn type to CDS if gene is selected
    ftype = 'CDS' if gene else ftype

    # Get a stream of
    recs = parser.get_records(fnames)

    if len(recs) < 1:
        utils.error("no sequence records found in data")

    recs = map(remapper, recs)

    # Keep all records when selecting by sequence id
    gff = not fasta
    if not seqid:
        feature_selection = (gene or type_ or protein or features or gff)
        if feature_selection:
            recs = filter(remove_source, recs)
        else:
            recs = filter(keep_source, recs)

    # Filter by seqid.
    recs = filter(seqid_selector(seqid), recs)

    # Filters gene and CDS
    recs = filter(gene_selector(gene), recs)

    # Extract proteins if requested.
    if protein:
        recs = filter(protein_filter, recs)
        recs = map(protein_extract, recs)

    # Apply additional filters.
    recs = filter(type_selector(ftype), recs)

    # Reverse complement the sequence
    recs = map(reverse_complement(revcomp), recs)

    # Slicing depends on output type.
    if fasta:
        recs = map(sequence_slicer(start=start, end=end), recs)
    else:
        recs = filter(interval_selector(start=start, end=end), recs)

    # Apply the translation
    recs = map(translate_recs(translate), recs)

    # Select the formatter.
    if fasta:
        # Fasta formatter.
        formatter = fasta_formatter
    else:
        # GFF formatter.
        print("##gff-version 3")
        formatter = gff_formatter

    # Display the results.
    for rec in recs:
        formatter(rec)
