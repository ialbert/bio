"""
Converts across formats
"""
import gzip
import os
import sys
from collections import OrderedDict, defaultdict
from itertools import count

from biorun import utils

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

# GenBank terms to remap according to sequence ontology.
SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}


def is_fasta(fname):
    fname = fname.lower()
    exts = "fa fasta fa.gz fasta.gz".split()
    found = filter(lambda x: fname.endswith(x), exts)
    return any(found)


def is_genbank(fname):
    fname = fname.lower()
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
        self.locations = [(loc.start, loc.end, loc.strand) for loc in
                          self.feat.location.parts] if feat is not None else []


def get_records(recs, format):
    """
    Returns sequence features
    """

    for obj in recs:

        # Handles Fasta input. TODO: This here needs to be rewritten
        if format == 'fasta':
            rec = SeqRecord(seq=obj.seq, name=obj.name, description=obj.description, id=obj.id)
            out = Record(rec=rec, feat=None, seqid=obj.id, annot={}, ftype=Record.SOURCE, strand=1, start=1,
                         end=len(obj.seq))
            yield out
            continue

        # Handles regular GenBank
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


def parse_stream(fname, format='genbank'):
    """
    Parses a filename with the appropriate readers.
    """

    if hasattr(fname, 'read'):
        stream = fname
    else:
        if not os.path.isfile(fname):
            logger.error(f"file not found: {fname}")
            sys.exit()
        stream = gzip.open(fname) if fname.endswith("gz") else open(fname)
        format = 'fasta' if is_fasta(fname) else format

    recs = SeqIO.parse(stream, format=format)

    recs = get_records(recs, format=format)

    return recs


def fasta_formatter(rec):
    print(rec.obj.format("fasta"), end='')


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


def reverse_complement(flag):
    """
    Reverse complement
    """

    def func(rec):
        if flag:
            rec.obj.seq = rec.obj.seq.reverse_complement()
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


def feature2gff(seqid, ftype, start, end, strand, uid, name, pid=None):
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
    data = [seqid, ".", ftype, start, end, ".", strand, phase, attr]

    return data


def gff_formatter(rec):
    """
    Formats a record as GFF.
    """
    # Parent feature
    data = feature2gff(start=rec.start, end=rec.end, ftype=rec.type, uid=rec.id, name=rec.name, strand=rec.strand,
                       seqid=rec.seqid, pid=None)
    line = "\t".join(map(str, data))

    # Parent id.
    pid = rec.id

    if rec.type == "mRNA":
        print(line)
        ftype = "exon"
    else:
        ftype = rec.type

    # Generate the locations
    for start, end, strand in rec.locations:
        name = rec.name
        uid = next(counter)

        data = feature2gff(start=start, end=end, ftype=ftype, uid=uid, name=name,
                           strand=strand, seqid=rec.seqid, pid=rec.id)
        line = "\t".join(map(str, data))
        print(line)


def run(features=False, protein=False, translate=False, fasta=False, revcomp=False,
        start='1', end=None, type_='', id_='', name='', gene='', alias=None, fnames=[]):
    """
    Converts data to different formats.
    """
    global ALIAS

    # Comes as tuple
    fnames = list(fnames)

    # Adding standard input as a stream
    if not sys.stdin.isatty():
        fnames.append(sys.stdin)

    # Generate the ALIAS remapping.
    ALIAS = utils.parse_alias(alias) if alias else {}

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

    # Selects sources only when no other feature specific option is set.
    source_flag = not (gene or name or type_ or translate or protein or features)

    # Use the source only in fasta mode.
    source_flag = source_flag and fasta

    # Select the formatter.
    if fasta:
        # Fasta formatter.
        formatter = fasta_formatter
    else:
        # GFF formatter.
        print("##gff-version 3")
        formatter = gff_formatter

    # Handle each input separately.
    for fname in fnames:

        # Parse the input into records.
        recs = parse_stream(fname)

        # Remap aliases.
        recs = map(remapper, recs)

        # Should we keep the source
        recs = filter(source_only(source_flag), recs)

        # Filter by sequence name
        recs = filter(name_selector(name), recs)

        # Filters gene and CDS
        recs = filter(gene_selector(gene), recs)

        # Filter by seqid.
        recs = filter(seqid_selector(seqid), recs)

        # Extract proteins if requested.
        if protein:
            recs = filter(protein_filter, recs)
            recs = map(protein_extract, recs)

        # Apply additional filters.
        recs = filter(type_selector(ftype), recs)

        if fasta:
            recs = map(sequence_slicer(start=start, end=end), recs)
        else:
            recs = filter(interval_selector(start=start, end=end), recs)

        # Reverse complement the sequence
        recs = map(reverse_complement(revcomp), recs)

        recs = map(translate_recs(translate), recs)

        # Display the results.
        for rec in recs:
            formatter(rec)
