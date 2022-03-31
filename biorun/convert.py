"""
Converts across formats
"""
import re, os
import sys
import string
from itertools import count
from biorun import utils, parser
import json

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
    fields = fields.split(",")
    hasN = 'N' in fields

    coll = ["{%s}" % f for f in fields]

    patt = "\t".join(coll)

    def func(rec):
        param = get_params(rec, hasN=hasN)
        try:
            print(patt.format(**param))
        except KeyError as exc:
            valid = ",".join(list(param.keys()))
            utils.error(f"Invalid formatter: {exc}", stop=False)
            utils.error(f"Valid values: {valid}")

    return func


def fasta_formatter(rec):
    seqrec = parser.make_seqrec(rec)
    print(seqrec.format("fasta"), end='')


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
        return rec.gene in genes if name else True

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
            cond = cpatt.search(rec.id) or cpatt.search(str(rec.desc))
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


COUNTER = count(1)


def get_params(rec, hasN=False):
    """
    Makes a dictionary out of parameters.
    """

    # text = json.dumps(ann, indent=4)
    # print (text)
    if hasN:
        nfunc = lambda rec: rec.seq.count('N')
    else:
        nfunc = lambda rec: 0

    params = dict(
        isolate=rec.get_parent_ann("isolate"),
        country=rec.get_parent_ann("country"),
        date=rec.get_parent_ann("collection_date"),
        pub_date=rec.get_parent_ann("date"),
        host=rec.get_parent_ann("host"),
        gene=rec.gene,
        type=rec.type,
        product=rec.product,
        size=len(rec.seq),
        source=rec.source,
        id=rec.id,
        count=next(COUNTER),
        N=nfunc(rec),
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

            # text = ascii(text)

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


def features_only(rec):
    return rec.type != parser.SOURCE


def protein_filter(rec):
    return "translation" in rec.ann


def protein_extract(rec):
    trans = parser.first(rec.ann, "translation")
    rec.seq = Seq(trans)
    return rec

def trim_maker(base):
    def trim_func(rec):
        rec.seq = rec.seq.rstrip(base)
        rec.seq = rec.seq.strip('N')
        return rec
    return trim_func
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


def feature2gff(source, ftype, start, end, strand, uid, name, pid=None):
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
    data = [source, ".", ftype, start, end, ".", strand, phase, attr]

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
                       name=rec.id, strand=rec.strand,
                       source=rec.source, pid=None)

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
        name = rec.id
        uid = next(parser.counter)

        data = feature2gff(start=start, end=end, ftype=ftype, uid=uid, name=name,
                           strand=strand, source=rec.source, pid=rec.id)
        line = "\t".join(map(str, data))
        print(line)


def run(protein=False, translate=False, fasta=False, revcomp=False, features=False, trim=False, olap='', table=False, fields='',
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

    # When to produce tabular output
    gff = not fasta

    # When are we producing features
    keepall = seqid or match or table

    features = features or gene or type_ or protein or gff

    # Keep the features or the sources only
    if not keepall:
        if features:
            recs = filter(features_only, recs)
        else:
            recs = filter(source_only, recs)

    if seqid:
        # Filter by seqid.
        recs = filter(seqid_selector(seqid), recs)

    if match:
        # Apply regular expression match.
        recs = filter(regex_selector(match), recs)

    if trim:
        # Get the trimming function.
        trim_func = trim_maker(trim)
        # Apply the trimming function.
        recs = map(trim_func, recs)

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
