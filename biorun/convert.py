"""
Converts across formats
"""
import re, os
import sys
import string

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


def table_formatter(fields):

    coll = ["{%s}" % f for f in fields.split(",")]

    patt = "\t".join(coll)

    def func(rec):
        param = get_params(rec)
        print(patt.format(**param))

    return func

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


def true(x):
    return True


def false(x):
    return True


def overlap_selector(olap):
    olaps = set(olap.split(","))

    if not olap:
        func = false
    else:
        try:
            values = list(map(int, olaps))
        except ValueError as exc:
            utils.error(f"Invalid numeric values: {olaps}")

        def func(rec):
            for value in values:
                if rec.start <= value <= rec.end:
                    return True
            return False

    return func


def type_selector(ftype):
    if ftype == 'all':
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


def regex_selector(patt):
    cpatt = re.compile(patt)

    def func(rec):
        if patt:
            cond = cpatt.search(rec.id) or cpatt.search(rec.description)
            return cond
        else:
            return True

    return func


def seqid_selector(seqid):
    targets = set(seqid.split(","))

    def func(rec):
        return rec.id in targets if seqid else True

    return func


def ascii(text):
    """
    Removes punctuation from text
    """
    remove = r"!\"$%&',/:;<>?@\^`{}~"
    table = str.maketrans('', '', remove)
    text = text.translate(table)
    return text

def get_params(rec):
    """
    Makes a dictionary out of parameters.
    """
    params = dict(
        isolate=rec.annot.get("isolate", ['.'])[0],
        country=rec.annot.get("country", ['.'])[0],
        date=rec.annot.get("collection_date", ['.'])[0],
        pub_date=rec.annot.get("date", '.'),
        host=rec.annot.get("host", ['.'])[0],
        gene=rec.gene or '.',
        type=rec.type,
        size=len(rec.seq),
        reference=rec.anchor,
        id=rec.id,
    )
    return params

def rename_sequence(patt):
    """
    Renames records base on a formatting pattern such as {isolate} or by a file.
    """

    if os.path.isfile(patt):
        # Renaming can be done from a file as well.
        mapper = utils.parse_alias(patt)

        def func(rec):
            rec.id = mapper.get(rec.id) or rec.id
            return rec
    else:

        # Allow controll characters.
        patt = bytes(patt, "utf-8").decode("unicode_escape")

        # Pattern based renames.
        def func(rec):

            params = get_params(rec)

            text = patt.format(**params)
            text = "_".join(text.split())

            #text = ascii(text)

            rec.id = text

            return rec

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


def source_only(rec):
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


def run(protein=False, translate=False, fasta=False, revcomp=False, source=False, olap='', table=False, fields='',
        start='1', end=None, type_='', id_='', match='', gene='', rename=None, fnames=[]):
    """
    Converts data to different formats.
    """

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

    # Keep the genome only
    if source:
        recs = filter(source_only, recs)

    if seqid:
        # Filter by seqid.
        recs = filter(seqid_selector(seqid), recs)

    if match:
        # Apply regular expression match.
        recs = filter(regex_selector(match), recs)

    # Select by overlap
    recs = filter(overlap_selector(olap), recs)

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

    # Apply a dynamic renaming
    if rename:
        recs = map(rename_sequence(rename), recs)

    # Select the formatter.
    if table:
        formatter = table_formatter(fields)
    elif fasta:
        # Fasta formatter.
        formatter = fasta_formatter
    else:
        # GFF formatter.
        print("##gff-version 3")
        formatter = gff_formatter

    # Display the results.
    for rec in recs:
        formatter(rec)
