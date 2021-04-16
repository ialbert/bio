"""
Represents a BioPython SeqRecords as JSON data structure.

There is functionality to convert from SeqRecord to JSON and vice versa.

Storing data as JSON instead of FASTA or GENBANK allows for much faster processing of the data.
"""
import sys, os, gzip, json
from collections import OrderedDict
from biorun import utils, const
from itertools import count
from collections import defaultdict
import plac
import itertools

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import Reference, CompoundLocation, FeatureLocation
except ImportError as exc:
    utils.error(f"{__name__}", stop=False)
    utils.error(f"{exc}", stop=False)
    utils.error(f"This software requires biopython.", stop=False)
    utils.error(f"Try: conda install biopython", stop=True)

logger = utils.logger

# Globally generates unqiue counters.
SEQ_COUNTER = count(1)

# Keeps track of unique items that belong to a name.
UNIQUE = defaultdict(int)

def first(item, key, default=""):
    """
    Shortcut to obtain the first element of the list
    """
    return item.get(key, [default])[0]


def json_ready(value):
    """
    Serializes elements in containers to a type that can be turned into JSON.
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

def get_next_count(label, ftype):
    """
    Counts
    """
    global UNIQUE
    key = f"{label}-{ftype}"
    UNIQUE[key] += 1
    return UNIQUE[key]

def fill_name(f):
    """
    Attempts to generate an unique id and a parent from a BioPython SeqRecord.
    Mutates the feature dictionary passed in as parameter.
    """
    global UNIQUE

    # Get the type
    ftype = f['type']

    # Get gene name
    gene_name = first(f, "gene")

    # Will attempt to fill in the uid from attributes.
    uid = ''

    # Deal with known types.
    if ftype == 'gene':
        name = gene_name or first(f, "locus_tag")
        uid = name
    elif ftype == 'CDS':
        count = get_next_count(ftype=ftype, label=gene_name)
        prot = first(f, "protein_id") or f"{gene_name}-CDS-{count}"
        uid = f"{prot}"
        name = prot
    elif ftype == 'mRNA':
        count = get_next_count(ftype=ftype, label=gene_name)
        uid = first(f, "transcript_id") or f"{gene_name}-mRNA-{count}"
        name = uid
    elif ftype == "exon":
        name = gene_name
    else:
        name = first(f, "organism") or first(f, "transcript_id") or None
        uid = first(f, "transcript_id")

    # Set the unique identifier.
    f['id'] = uid or f"{ftype}-{next(SEQ_COUNTER)}"

    # Set the feature name.
    f['name'] = name or ftype

    return f



def make_attr(feat, color=None):
    """
    Creates GFF style attributes from JSON fields.
    """

    # The feature name.
    name = feat['name']
    uid = feat['id']
    pid = feat.get("parent_id")

    data = [f"ID={uid}"]

    # Add parent id
    if pid:
        data.append(f"Parent={pid}")

    # Add the data name.
    data.append(f"Name={name}")

    # Other gff attributes.
    pairs = [(k, v) for k, v in feat.items() if k not in const.SKIP_GFF_ATTR]

    # Fill in all the fields.
    for key, value in pairs:
        data.append(f"{key}={value[0]}")

    # Attach a color to the feature.
    if color:
        data.append(f"color={color}")

    return ";".join(data)


def rec_desc(f):
    """
    Generates a record description from JSON feature.
    """
    return make_attr(f)

def convert_genbank(recs, seqid=None):
    """
    Converts BioPython SeqRecords obtained from a GENBANK file to a JSON representation.
    """
    global SEQ_COUNTER

    # The outer dictionary containing multiple records.
    data = []

    count = 0
    # Add each record separately.
    for rec in recs:
        count += 1
        # Each individual SeqRecords is a dictionary.
        item = dict()

        # Fill the standard SeqRecord fields.
        item[const.SEQID] = seqid or rec.id
        item[const.DEFINITION] = rec.description
        item[const.DBLINK] = rec.dbxrefs
        item[const.LOCUS] = rec.name
        item[const.FEATURE_COUNT] = len(rec.features)
        item[const.ORIGIN_SIZE] = len(rec.seq)

        # Fill in all annotations.
        for key, value in rec.annotations.items():
            item[key] = json_ready(value)

        # Fill in the features.
        feats = []
        for feat in rec.features:

            # Feature type.
            ftype = feat.type

            # Remap GenBank terms to Sequence Ontology terms.
            ftype = const.SEQUENCE_ONTOLOGY.get(ftype, ftype)

            # Feature strand.
            strand = feat.strand

            # Feature coordinates are 1 based.
            start = int(feat.location.start) + 1
            end = int(feat.location.end)

            # Location operator.
            oper = feat.location_operator

            # Location coordinates are 1 based.
            location = [(loc.start + 1, loc.end, loc.strand) for loc in feat.location.parts]

            # Feature attributes
            attrs = dict(type=ftype, name='', id='', start=start, end=end, strand=strand, location=location,
                         operator=oper)

            # Fills in the additional qualifiers.
            for (k, v) in feat.qualifiers.items():
                attrs[k] = json_ready(v)

            # Adjust uid, parent and name using the attributes.
            attrs = fill_name(attrs)

            # Append the attributes as a record.
            feats.append(attrs)

        # Add the features
        item[const.FEATURES] = feats

        # Save the sequence as well
        item[const.ORIGIN] = str(rec.seq)

        # Add the item to the list.
        data.append(item)

    return data


def convert_fasta(recs, seqid=None):
    """
    Converts BioPython SeqRecords obtained from a FASTA file to JSON representation.
    """
    data = []
    for rec in recs:
        item = dict()
        item[const.SEQID] = seqid or rec.id
        item[const.LOCUS] = rec.name
        item[const.DEFINITION] = rec.description
        item[const.ORIGIN] = str(rec.seq)
        item[const.FEATURES] = []
        data.append(item)
    return data

def make_jsonrec(seq, seqid=None):
    """
    Makes a simple JSON representation for a text
    """
    global SEQ_COUNTER
    name = f"{next(SEQ_COUNTER)}"
    uid = seqid or f"{name}"
    data = []
    item = dict()
    item[const.SEQID] = uid
    item[const.LOCUS] = ''
    item[const.DEFINITION] = ''
    item[const.ORIGIN] = str(seq)
    start, end, strand = 1, len(seq), 1
    oper = None,
    location = [[start, end, strand]]
    ftype = "sequence"
    attrs = dict(start=start, end=end, type=ftype, strand=strand, location=location, operator=oper,
                 name=uid, id=uid)
    item[const.FEATURES] = [
        attrs
    ]
    data.append(item)

    return data


def origin_records(item):
    """
    Returns the origin sequence from an JSON item
    """
    # Get origin sequence.
    text = item[const.ORIGIN]

    # Transform to BioPython Sequence.
    seq = Seq(text)

    # Additional sequence attributes.
    seqid = item[const.SEQID]
    locus = item[const.LOCUS]
    desc = item.get("definition", "")

    # Build a BioPython sequence record.
    rec = SeqRecord(seq, id=seqid, name=locus, description=desc)

    # Must be iterable
    yield rec


def feature_records(data, translate=False):
    """
    Yields BioPython SeqRecords from JSON data.
    """

    # The feature generator
    feats = data[const.FEATURES]

    # We can extract DNA sequences from this if needed.
    origin = data[const.ORIGIN]

    # Create a SeqRecord from each feature.
    for f in feats[1:]:

        ftype = f['type']

        # Skip the origin.
        if ftype == 'region':
            continue

        # Concatenate locations
        locations = f.get("location", [])

        # Build the sequence for the record.
        dna = Seq('')
        for x, y, strand in locations:
            chunk = Seq(origin[x - 1:y])
            if strand == -1:
                chunk = chunk.reverse_complement()
            dna += chunk

        # Generate a simple description.
        desc = simple_description(f)

        # Build the sequence record.
        rec = SeqRecord(dna, id=f['id'], description=desc)

        yield rec


def simple_description(f):
    return f"{f['type']}"

def protein_records(item):
    """
    Yields BioPython SeqRecords for JSON features that do have protein translations.
    """


    # All features for the item
    feats = item[const.FEATURES]

    # Filtering function for translations.
    has_translation = lambda f: f.get('translation', [''])[0]

    # Features with translation.
    feats = filter(has_translation, feats)

    # Produce the translation records.
    for f in feats:
        # Fetch the translation.
        trans = first(f, "translation")

        # Create BioPython sequence.
        seq = Seq(trans)

        # Generate a simple description
        desc = simple_description(f)

        # Form the sequence record.
        rec = SeqRecord(seq, id=f['id'], description=desc)

        yield rec


def parse_stream(stream, type='genbank'):
    """
    Parses a recognized file into a JSON representation
    """

    # Cascade over the known formats.
    if type=='fasta':
        recs = SeqIO.parse(stream, format=const.FASTA)
        data = convert_fasta(recs)
    else:
        recs = SeqIO.parse(stream, format=const.GENBANK)
        data = convert_genbank(recs)

    text = json.dumps(data, indent=4)
    print(text)

@plac.opt("type_", "input type", choices=["gb", "fasta"])
def run(type_="gb"):
    """
    The main JSON converter.
    """
    parse_stream(sys.stdin)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
