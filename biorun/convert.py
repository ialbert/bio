"""
Converts across formats
"""
import gzip
import os
import sys
from collections import OrderedDict, defaultdict
from itertools import count

import biorun.libs.placlib as plac
from biorun import utils
from biorun.alias import ALIAS

# Module level logger.
logger = utils.logger

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import Reference, CompoundLocation, FeatureLocation
except ImportError as exc:
    logger.error(f"{__name__}", stop=False)
    logger.error(f"{exc}", stop=False)
    logger.error(f"This software requires biopython.")
    logger.error(f"Try: conda install biopython>=1.78")
    sys.exit()

# GenBank terms to remap according to sequence ontology.
SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}


def is_fasta(fname):
    exts = "fa fasta fa.gz fasta.gz".split()
    found = filter(lambda x: fname.endswith(x), exts)
    return any(found)


def is_genbank(fname):
    exts = "gb gb.gz genbank genbank.gz gpff gpff.gz".split()
    found = filter(lambda x: fname.endswith(x), exts)
    return any(found)


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


def first(data, key, default=""):
    # First element of a list value that is stored in a dictionary by a key.
    return data.get(key, [default])[0]


counter = count(1)

UNIQUE = defaultdict(int)


def next_count(ftype):
    UNIQUE[ftype] += 1
    return f'{ftype}-{UNIQUE[ftype]}'


def guess_name(ftype, annot):
    """
    Attempts to generate an unique id name for a BioPython feature
    """
    uid = desc = ''

    if ftype == 'gene':
        name = first(annot, "gene")
        desc = first(annot, "locus_tag")
    elif ftype == 'CDS':
        name = first(annot, "protein_id")
        desc = first(annot, "product")
    elif ftype == 'mRNA':
        name = first(annot, "transcript_id")
        desc = ftype
    elif ftype == "exon":
        name = first(annot, "gene")
        uid = next_count(ftype)
    else:
        name = next_count(ftype)
        desc = first(annot, "product")

    desc = f"{ftype} {desc}"
    name = name or next_count(f"unknown-{ftype}")
    uid = uid or name
    return uid, name, desc


class RecAttrs:
    """
    Attributes that describe a grouping of records.
    """

    def __init__(self, obj, **kwds):
        self.obj = obj
        self.id = obj.id
        self.name = obj.name
        self.seq = obj.seq
        self.title = obj.description
        pairs = [(k, json_ready(v)) for (k, v) in self.obj.annotations.items()]
        self.annot = dict(pairs)


class Record:
    """
    Unified representation of BioPython SeqRecord
    """
    SOURCE = "source"

    def __init__(self, rec, seqid, feat, ftype, start, end, strand, annot={}, attrs=None):
        # The original object
        self.obj = rec
        self.seqid = ALIAS.get(seqid, seqid)
        self.feat = feat
        self.attrs = attrs
        self.type = SEQUENCE_ONTOLOGY.get(ftype, ftype)
        self.id = rec.id
        self.name = rec.name
        self.title = rec.description

        self.start, self.end, self.strand = start, end, strand

        # Fill the annotations
        self.annot = annot

        # Store the locations.
        self.locations = [(loc.start, loc.end, loc.strand) for loc in self.feat.location.parts]


def get_records(recs):
    """
    Returns sequence features
    """
    for obj in recs:

        for feat in obj.features:
            # Normalize the feature type.
            ftype = SEQUENCE_ONTOLOGY.get(feat.type, feat.type)

            # Sequence for this feature
            seq = feat.extract(obj.seq)

            # The start/end locations
            start, end = int(feat.location.start), int(feat.location.end)

            # Qualifiers are transformed into annotations.
            pairs = [(k, json_ready(v)) for (k, v) in feat.qualifiers.items()]
            annot = dict(pairs)

            # Create an id, name and descriptions
            uid, name, desc = guess_name(ftype=ftype, annot=annot)

            # Source objects need to be matched to top level annotations
            if ftype == Record.SOURCE:
                uid = obj.id
                name = obj.name
                desc = obj.description
                pairs = [(k, json_ready(v)) for (k, v) in obj.annotations.items()]
                annot.update(pairs)

            # Remap the uid
            uid = ALIAS.get(uid, uid)

            # Create the sequence record.
            rec = SeqRecord(seq=seq, name=name, description=desc, id=uid)

            # A Wrapper object that mixes Seqrecord and Features
            out = Record(rec=rec, feat=feat, seqid=obj.id, annot=annot, ftype=ftype, strand=feat.strand, start=start,
                         end=end)

            yield out


def parse(fname):
    """
    Parses a filename with the appropriate readers.
    """
    stream = gzip.open(fname) if fname.endswith("gz") else open(fname)
    if is_fasta(fname):
        recs = SeqIO.parse(stream, format="fasta")
        recs = get_records(recs)
    elif is_genbank(fname):
        recs = SeqIO.parse(stream, format="genbank")
        recs = get_records(recs)
    else:
        logger.error(f"file extension not recognized: {fname}")
        sys.exit()

    return recs


def make_record(text, seqid=1, locus="", desc=""):
    seq = Seq(text)
    recs = [SeqRecord(seq=seq, id=seqid, name=locus, description=desc)]
    recs = get_records(recs)
    return recs


def read_input(fname, store=None, interactive=False):
    """
    Attempts load the correct input.
    """

    # Item is a valid file.
    if os.path.isfile(fname):
        recs = parse(fname)
        return recs

    # Generate an interactive data.
    if interactive:
        recs = make_record(text=fname)
        return recs

    # Invalid data.
    logger.error(f"file not found: {fname}")
    sys.exit()


def fasta_formatter(rec):
    print(rec.obj.format("fasta"), end='')


from biorun.gff import gff_formatter


def remapper(rec):
    rec.id = ALIAS.get(rec.id, rec.id)
    return rec


def sequence_slicer(start=0, end=None):
    def func(rec):
        if start or end:
            seqlen = len(rec.obj.seq)

            endx = seqlen if end is None else end
            endx = seqlen + endx if endx < 0 else endx

            # Zero based shift
            if start < 0:
                startx = seqlen + start + 1
            else:
                startx = start + 1

            rec.obj.description = f"{rec.title} [{startx}:{endx}]"
            rec.obj.seq = rec.obj.seq[start:end]
        return rec

    return func


def type_selector(ftype):
    types = set(ftype.split(","))

    def func(rec):
        return rec.type in types if ftype else True

    return func


def gene_selector(name):
    genes = set(name.split(","))

    def func(rec):
        return first(rec.annot, "gene") in genes if name else True

    return func


def name_selector(name):
    names = set(name.split(","))

    def func(rec):
        return rec.id in names if name else True

    return func


def seqid_selector(seqid):
    targets = set(seqid.split(","))

    def func(rec):
        return rec.seqid in targets if seqid else True

    return func


def translate_recs(flag):
    def func(rec):
        if flag:
            endx = len(rec.obj.seq) // 3
            rec.obj.seq = rec.obj.seq[0:endx * 3].translate()
        return rec

    return func


def source_only(flag):

    def func(rec):
        return rec.type == Record.SOURCE if flag else True

    return func


def protein_filter(rec):
    return "translation" in rec.annot


def protein_extract(rec):
    rec.obj.seq = Seq(rec.annot.get("translation")[0])
    return rec


@plac.pos("data", "input data")
@plac.flg("features", "convert the features", abbrev='F')
@plac.flg("fasta_", "convert to fasta")
@plac.flg("gff_", "convert to gff")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "filter for a sequence id")
@plac.opt("name", "filter for a sequence name")
@plac.opt("gene", "filter for a gene name", abbrev='G')
@plac.flg("protein", "operate on the protein sequences", abbrev='P')
@plac.flg("translate", "translate DNA sequences", abbrev='R')
def run(features=False, protein=False, translate=False, gff_=False, fasta_=False,
        start='1', end=None, type_='', id_='', name='', gene='', *fnames):
    """
    Convert data to various formats
    """

    # When to generate features
    features = features or (type_ or gene or translate or name)

    # Parse start and end into user friendly numbers.
    start = utils.parse_number(start)

    # Positive coordinates moved to 0 based.
    if start >= 0:
        start = start - 1

    end = utils.parse_number(end)
    ftype = type_
    seqid = id_

    # Turn type to CDS if gene is selected
    ftype = 'CDS' if gene else ftype

    # Default format is fasta if nothing is specified.
    fasta = False if (gff_ and not fasta_) else True

    # Selects sources only when no other feature specific option is set.
    source_flag = not(gene or name or type_ or translate or protein)

    # GFF mode produces all features
    source_flag = False if gff_ else source_flag

    # Select the formatter.
    if fasta:
        formatter = fasta_formatter
    else:
        print("##gff-version 3")
        formatter = gff_formatter

    # Handle each input separately.
    for fname in fnames:

        # Produces the input as a record generator.
        recs = read_input(fname, interactive=False)

        # Remap aliases/
        recs = map(remapper, recs)

        # Filter by sequence name
        recs = filter(name_selector(name), recs)

        # Should we keep the source
        recs = filter(source_only(source_flag), recs)

        # Filters gene and CDS
        recs = filter(gene_selector(gene), recs)

        # Filter by seqid
        recs = filter(seqid_selector(seqid), recs)

        # Extract proteins.
        if protein:
            recs = filter(protein_filter, recs)
            recs = map(protein_extract, recs)

        # Apply additonal filters.
        recs = filter(type_selector(ftype), recs)
        recs = map(sequence_slicer(start=start, end=end), recs)
        recs = map(translate_recs(translate), recs)

        # Display the results.
        for rec in recs:
            formatter(rec)

    sys.exit()
