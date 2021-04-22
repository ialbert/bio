"""
Represents a BioPython SeqRecords as JSON data structure.

There is functionality to convert from SeqRecord to JSON and back.

Storing data as JSON instead of FASTA or GENBANK allows for much faster processing.
"""
import sys, os, json, gzip, pathlib
from collections import OrderedDict
from collections import defaultdict
from functools import partial
from itertools import count
from biorun.alias import ALIAS

import plac
from biorun import utils, reclib

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

# JSON specific constants
DEFINITION = "definition"
DBLINK = "dblink"
LOCUS = "locus"
ID = "id"
ACCESSION = "accession"

# New keys not originally present in the GenBank record.
ORIGIN = "origin"
FEATURE_LIST = "features"
FEATURE_COUNT = "feature_count"
ORIGIN_LEN = "sequence_len"
TAXIDS = "taxids"

# GenBank terms to remap according to sequence ontology.
SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}

# Globally generates unqiue counters.
SEQ_COUNTER = count(1)

# Keeps track of unique items that belong to a name.
UNIQUE = defaultdict(int)


def find_taxid(rec):
    """
    Attempts to extract the taxonomy id from a record.

    For GenBank files the taxid is the listed in the first feature (the source).
    In the JSON file it will be represented like so.

    "db_xref": [
        "taxon:11191"
    ],

    """
    try:
        feats = rec.get(FEATURE_LIST, [])
        db_xref = feats[0].get("db_xref", []) if feats else []
        taxids = [x for x in db_xref if x.startswith("taxon:")]
    except Exception as exc:
        logger.error(f"Error parsing taxid: {exc}")
        taxids = []

    return taxids


def first(item, key, default=""):
    """
    First element of a list value, stored in a dictionary by a key.
    """
    return item.get(key, [default])[0]


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


def get_next_count(name, ftype):
    """
    Keep individual counts for elements named certain way.
    """
    global UNIQUE
    key = f"{name}-{ftype}"
    UNIQUE[key] += 1
    return UNIQUE[key]


def fill_name(f, seqid=''):
    """
    Attempts to generate an unique id and a parent from a BioPython SeqRecord.
    Mutates the feature dictionary passed in as parameter and adds a
    """

    # Get the type
    ftype = f['type']

    # Get gene name
    gene_name = first(f, "gene")

    uid = name = ''

    # Deal with known types.
    if ftype == 'source':
        uid = seqid
        name = seqid.split(".")[0]
    elif ftype == 'gene':
        name = gene_name or first(f, "locus_tag")
        uid = name
    elif ftype == 'CDS':
        count = get_next_count(ftype=ftype, name=gene_name)
        prot = first(f, "protein_id") or f"{gene_name}-CDS-{count}"
        uid = f"{prot}"
        name = prot
    elif ftype == 'mRNA':
        count = get_next_count(ftype=ftype, name=gene_name)
        uid = first(f, "transcript_id") or f"{gene_name}-mRNA-{count}"
        name = uid
    elif ftype == "exon":
        name = gene_name
    else:
        name = first(f, "organism") or first(f, "transcript_id") or None
        uid = first(f, "transcript_id")

    # Make a copy of old feature.
    feat = dict(f)

    # Set the unique identifier.
    feat['id'] = uid or f"{ftype}-{next(SEQ_COUNTER)}"

    # Set the feature name, fallback to feature type
    feat['name'] = name or ftype

    return feat


def convert_genbank(recs, remap={}):
    """
    Converts BioPython SeqRecords obtained from a GENBANK file to a JSON representation.
    """
    global SEQ_COUNTER

    # The return data will be a list of records.
    data = []

    # Add each record into the json data.
    for rec in recs:

        # Each SeqRecords will be a dictionary.
        item = dict()

        # Fill the standard SeqRecord fields.
        item[ID] = ALIAS.get(rec.id, rec.id)
        item[ACCESSION] = rec.id
        item[DEFINITION] = rec.description
        item[DBLINK] = rec.dbxrefs
        item[LOCUS] = rec.name
        item[FEATURE_COUNT] = len(rec.features)
        item[ORIGIN_LEN] = len(rec.seq)
        item[TAXIDS] = []
        # Fill in all annotations.
        for key, value in rec.annotations.items():
            item[key] = json_ready(value)

        # Fill in the features.
        feats = []
        for feat in rec.features:

            # Feature type.
            ftype = feat.type

            # Remap GenBank terms to Sequence Ontology terms.
            ftype = SEQUENCE_ONTOLOGY.get(ftype, ftype)

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

            feat = fill_name(attrs, seqid=item[ID])

            # Append the attributes as a record.
            feats.append(feat)

        # Add the features
        item[FEATURE_LIST] = feats

        # Parse the taxonomic id
        item[TAXIDS] = find_taxid(item)

        # Save the sequence as well
        item[ORIGIN] = str(rec.seq)

        # Add the item to the list.
        data.append(item)

    return data


def convert_fasta(recs, remap={}):
    """
    Converts BioPython SeqRecords obtained from a FASTA file to JSON representation.
    """
    data = []
    for rec in recs:
        item = dict()
        item[ID] = remap.get(rec.name) or rec.id
        item[LOCUS] = rec.name
        item[DEFINITION] = rec.description
        item[ORIGIN] = str(rec.seq)
        item[FEATURE_LIST] = []
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
    item[ID] = uid
    item[LOCUS] = ''
    item[DEFINITION] = ''
    item[ORIGIN] = str(seq)
    start, end, strand = 1, len(seq), 1
    oper = None,
    location = [[start, end, strand]]
    ftype = "sequence"
    attrs = dict(start=start, end=end, type=ftype, strand=strand, location=location, operator=oper,
                 name=uid, id=uid)
    item[FEATURE_LIST] = [
        attrs
    ]
    data.append(item)

    return data


def origin_records(item):
    """
    Returns the origin sequence from an JSON item
    """
    # Get origin sequence.
    text = item[ORIGIN]

    # Transform to BioPython Sequence.
    seq = Seq(text)

    # Additional sequence attributes.
    seqid = item[ID]
    locus = item[LOCUS]
    desc = item.get("definition", "")

    # Build a BioPython sequence record.
    rec = SeqRecord(seq, id=seqid, name=locus, description=desc)

    # We will add a new attribute here.
    rec.__biotype__ = ORIGIN

    # Must be iterable
    yield rec


def feature_records(data, skip=True):
    """
    Yields BioPython SeqRecords from JSON data.
    """

    # The feature generator
    feats = data[FEATURE_LIST]

    # We can extract DNA sequences from this if needed.
    origin = data[ORIGIN]

    # Create a SeqRecord from each feature.
    for f in feats:

        # Type of the feature.
        ftype = f['type']

        # Skip the origin.
        if skip and ftype == 'source':
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
        rec = SeqRecord(dna, id=f['id'], description=desc, name=f['name'])

        # We will add a new attribute here.
        rec.__biotype__ = f['type']
        rec.__gene__ = first(f, "gene")

        yield rec


def simple_description(f):
    desc = f"type={f['type']}"
    gene = first(f, "gene")
    if gene:
        desc = f"{desc} gene={gene}"
    return desc


def protein_records(item):
    """
    Yields BioPython SeqRecords for JSON features that do have protein translations.
    """

    # All features for the item
    feats = item[FEATURE_LIST]

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
        rec = SeqRecord(seq, id=f['id'], description=desc, name=f['name'])

        # We will add a new attribute here.
        rec.__biotype__ = f['type']

        yield rec


def parse_stream(stream, type='genbank', remap={}):
    """
    Parses a recognized file into a JSON representation
    """

    # Cascade over the known formats.
    if type == "fasta":
        recs = SeqIO.parse(stream, format="fasta")
        data = convert_fasta(recs, remap=remap)
    else:
        recs = SeqIO.parse(stream, format="genbank")
        data = convert_genbank(recs, remap=remap)

    return data


def select_features(data, start=0, end=0, ftype=None, seqid=None, name=None, gene=None, match=None):
    """
    Selects features with certain properties
    """
    for item in data:

        itemid = item['id']

        # parent_id = gene_id = transcript_id = None
        features = item[FEATURE_LIST]

        if ftype:
            features = filter(lambda x: x['type'] == ftype, features)

        if seqid:
            if itemid != seqid:
                continue

        for feat in features:

            feat_id = feat['id']
            feat_type = feat['type']

            locations = feat["location"]
            size = len(locations)
            # Deal with multilevel locations
            if feat_type == "mRNA" or feat_type == "CDS":

                transcript_id = feat.get("transcript_id") or [feat_id]
                gene_id = feat.get("locus_tag") or feat.get("gene") or [feat_id]

                if feat_type == "mRNA":
                    child_type = "exon"
                    yield itemid, feat
                else:
                    child_type = "CDS"

                counter = count(1)

                for start, end, strand in locations:
                    uid = f"{child_type}-{feat_id}-{next(counter)}"
                    entry = dict(parent_id=feat_id, type=child_type, strand=strand,
                                 start=start, end=end, name=uid, id=uid,
                                 transcript_id=transcript_id, gene_id=gene_id)

                    yield itemid, entry
            elif size > 1:
                for start, end, strand in locations:
                    entry = dict(feat)
                    entry['parent_id'] = feat_id
                    entry['start'], entry['end'], entry['strand'] = start, end, strand
                    yield itemid, entry
            else:
                yield itemid, feat


def select_records(data, features=False, proteins=False,
                   translate=False, start=0, end=0, ftype=None, seqid=None, name=None, gene=None, match=None):
    """
    Selects records with certain properties
    """

    for item in data:

        # Turn on features when selecting for types
        features = True if (ftype or seqid or name or gene) else features

        # When to skip origin FASTA records.
        skip = False if (seqid or name or gene) else True

        if proteins:
            recs = protein_records(item)
        elif features:
            recs = feature_records(item, skip=skip)
        else:
            recs = origin_records(item)

        # Filter for a name
        func_gene = partial(reclib.filter_gene, gene=gene)
        recs = filter(func_gene, recs)

        # Filter for a name
        func_name = partial(reclib.filter_name, name=name)
        recs = filter(func_name, recs)

        # Filter for a sequence id
        func_id = partial(reclib.filter_seqid, seqid=seqid)
        recs = filter(func_id, recs)

        # Filter for a type
        func_type = partial(reclib.filter_type, ftype=ftype)
        recs = filter(func_type, recs)

        # Apply the slices.
        func_slice = partial(reclib.slice, start=start, end=end)
        recs = map(func_slice, recs)

        # Apply the translation.
        if translate:
            recs = map(reclib.translate, recs)

        # Yield all records in the data
        for rec in recs:
            yield rec


@plac.opt("alias", "alias file to rename sequences")
@plac.opt("type_", "input type", choices=["gb", "fasta"])
def run(type_="gb", alias=''):
    """
    The main JSON converter.
    """

    remap = utils.alias_dict(fname=alias)
    data = parse_stream(sys.stdin, remap=remap)
    text = json.dumps(data, indent=4)
    print(text)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
