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

from pprint import pprint
import itertools

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import Reference, CompoundLocation, FeatureLocation
except ImportError as exc:
    print(f"*** {exc}", file=sys.stderr)
    print(f"*** This software requires biopython.", file=sys.stderr)
    print(f"*** Try: conda install biopython", file=sys.stderr)
    sys.exit(1)

logger = utils.logger

COUNTER = count(1)
SEQNAME = count(1)

# Using this to keep track of unique items that belong to a name.
UNIQUE = defaultdict(int)


# Allow for resetting the global counter. Needed for keeping test labeled consistenly.
def reset_counter():
    global COUNTER, UNIQUE
    COUNTER = count(1)
    UNIQUE = defaultdict(int)


def reset_sequence_names():
    global SEQNAME
    reset_counter()

    # Adds two default sequences.
    SEQNAME = itertools.chain(["TARGET", "QUERY"], count(1))


def has_feature(item, name="gene"):
    """
    A filtering function to checks if a record contains keys with a name.
    """
    return item.get(name, [''])[0]


def filter_features(items, param, droporigin=False):
    """
    Filters features based on various parameters.
    """

    # Remove source as a valid feature.
    if droporigin:
        items = filter(lambda f: f.get('type') != 'region', items)

    # Filter by type.
    if param.type and param.type != "all":
        valid = set(map(lambda x: x.lower(), param.type.split(",")))
        valid = set(valid)
        items = filter(lambda f: f.get('type', '').lower() in valid, items)

    # Filter by gene.
    if param.gene:
        items = filter(lambda f: param.gene in f.get("gene", []), items)

    # Filter by name.
    if param.name:
        items = filter(lambda f: param.name in f.get("name", []), items)

    # Filter by external attributes.
    if param.match_field:
        items = filter(lambda f: param.match_value in f.get(param.match_field, []), items)

    # Filter by coordinates. Fasta files are cut by sequence.
    if (param.start or param.end) and not param.fasta:
        param.end = sys.maxsize if param.end is None else param.end
        items = filter(lambda f: param.start <= f.get('end') and param.end >= f.get('start'), items)

    # Filters by matching a regular expression
    if param.regexp:
        items = filter(lambda f: param.regexp.search(str(f)), items)

    return items


def find_taxid(rec):
    """
    Attempts to extract the taxonomy id from a record.

    For GenBank files the taxid is the listed in the first feature (the source)

    "db_xref": [
                "taxon:11191"
    ],

    """
    feats = rec.get(const.FEATURES, [])
    db_xref = feats[0].get("db_xref", []) if feats else []
    values = [x for x in db_xref if x.startswith("taxon:")]
    taxids = [x.split(":")[1] for x in values]
    return taxids


def first(item, key, default=""):
    """
    Shortcut to obtain the first element of the list
    """
    return item.get(key, [default])[0]


def rec_name(f):
    """
    Generates a record name from a JSON feature.
    """

    name = f['name']
    return name


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


def get_translation_records(item, param):
    """
    Yields BioPython SeqRecords for JSON features that do have translations.
    """

    # All features for the item
    feats = item[const.FEATURES]

    # Filtering function for translations.
    has_translation = lambda f: f.get('translation', [''])[0]

    # Features with translation.
    feats = filter(has_translation, feats)

    # Additional filters that may have been passed.
    feats = filter_features(feats, param=param)

    # Hoist the variables out.
    start, end = param.start, param.end

    # Produce the translation records.
    for f in feats:
        # Fetch the translation.
        trans = first(f, "translation")[start:end]

        # Set up the metadata.
        name = param.seqid or rec_name(f)

        # Create the sequence descriptor.
        desc = rec_desc(f)

        # Generate sequence record.
        seq = Seq(trans)

        # Form the sequence record.
        rec = SeqRecord(seq, id=name, description=desc)

        yield rec


def get_json_features(data):
    """
    Generates features from data.
    It will also generate the parent child relationships for hierachical features (mRNA and CDS).
    """

    feats = data[const.FEATURES]

    counter = count(1)

    for feat in feats:

        # The type of the location
        ftype = feat['type']

        # The unique id of the feature.
        feature_id = feat['id']

        # Overriden during hierachical data.
        location_type = ftype

        # Used for hierachical data.
        parent_id = None

        gene_id = transcript_id = None

        # Hierarchical data produce a parent feature.
        if ftype == 'mRNA' or ftype == 'CDS':

            # Child nodes will track parent.
            parent_id = feature_id

            # The GenBank nomenclature does not follow the Sequence Ontology,
            # CDS and exons are both considered to be part of mRNA.
            # In our model we reparent exons to transcripts and CDS to mRNA regions.

            # Set the parent/child types.
            if ftype == 'mRNA':
                # Keep the mRNA as parent of exons.
                parent_type, location_type = 'transcript', 'exon'
            else:
                # Make a new parent for CDS regions.
                parent_type, location_type = 'mRNA', 'CDS'

            # Figure out the gene and transcripts ids.
            gene_id = feat.get("locus_tag") or feat.get("gene") or parent_id
            transcript_id = feat.get("protein_id") or parent_id

            # Copy the attributes
            parent_feat = dict(feat)
            parent_feat['type'] = parent_type

            yield parent_feat

        # Generate JSON record for each location separately.
        for start, end, strand in feat["location"]:
            loc_feat = dict(feat)
            if parent_id:
                # Needs a new ID as the parent keeps the original id.
                loc_feat["id"] = f"{location_type}-{next(counter)}"

                # Assign the parent.
                loc_feat['parent_id'] = parent_id

                # Add a gene and transcript id attributes.
                if gene_id:
                    loc_feat['gene_id'] = gene_id
                if transcript_id:
                    loc_feat['transcript_id'] = transcript_id

            loc_feat['type'] = location_type
            loc_feat['start'], loc_feat['end'], loc_feat['strand'] = start, end, strand
            yield loc_feat


def get_feature_records(data, param):
    """
    Yields BioPython SeqRecords from JSON data.
    """

    # The feature generator
    feats = data[const.FEATURES]

    # Filter the features.
    feats = filter_features(feats, param=param, droporigin=True)

    # We can extract DNA sequences from this if needed.
    origin = data[const.ORIGIN]

    # Shortcuts to coordinates.
    start, end = param.start, param.end

    # Create a SeqRecord from each feature.
    for f in feats:

        # Make a name
        name = rec_name(f)

        # Concatenate locations
        locations = f.get("location", [])

        # Build the sequence for the record.
        dna = Seq('')
        for x, y, strand in locations:
            chunk = Seq(origin[x - 1:y])
            if strand == -1:
                chunk = chunk.reverse_complement()
            dna += chunk

        # Initialize the description
        desc = [f['type']]

        # Figure out description for slices.
        if start or end:
            _end = len(dna) if end is None else end
            desc.append(f'[{start + 1}:{_end}]')

        # Slice the resulting DNA sequence.
        seq = dna[start:end]

        try:
            # Perform sequence transformation if needed.
            if param.revcomp:
                seq = seq.reverse_complement()
                desc.append("reverse complemented")

            if param.reverse:
                seq = seq[::-1]
                desc.append("reversed")

            if param.complement:
                seq = seq.complement()
                desc.append("complemented")

            if param.translate:
                seq = seq.translate()
                desc.append("translated")

            if param.transcribe:
                seq = seq.transcribe()
                desc.append("transcribed DNA")

        except Exception as exc:
            utils.error(exc)

        # Make a description
        desc = " ".join(desc) if desc else rec_desc(f)

        # Build the sequence record.
        rec = SeqRecord(seq, id=name, description=desc)

        # Sanity check for translation
        if param.translate:
            expected = first(f, "translation")
            observed = str(seq)
            # Stop codon is present in the CDS but not in the translation.
            if expected and expected != observed[:-1]:
                logger.info(f"translation mismatch for: {rec.id}")

        yield rec


def get_origin(item, param):
    """
    Returns the origin sequence from an JSON item
    """
    # Prints the source sequence
    text = item[const.ORIGIN][param.start:param.end]
    seq = Seq(text)

    # Translates the origin.
    if param.translate:
        seq = seq.translate() if param.translate else seq

    # Fill the sequence attributes.
    desc = item[const.DEFINITION]
    seqid = item[const.SEQID]
    locus = item[const.LOCUS]
    seqid = param.seqid or seqid

    # Create the sequence record.
    rec = SeqRecord(seq, id=seqid, name=locus, description=desc)

    yield rec


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
        name = gene_name
        uid = first(f, "locus_tag") or gene_name
    elif ftype == 'CDS':
        count = get_next_count(ftype=ftype, label=gene_name)
        uid = first(f, "protein_id") or f"{gene_name}-CDS-{count}"
        name = uid
    elif ftype == 'mRNA':
        count = get_next_count(ftype=ftype, label=gene_name)
        uid = first(f, "transcript_id") or f"{gene_name}-mRNA-{count}"
        name = uid
    elif ftype == "exon":
        name = gene_name
    else:
        name = first(f, "organism") or None

    # Set the unique identifier.
    f['id'] = uid or f"{ftype}-{next(COUNTER)}"

    # Set the feature name.
    f['name'] = name or ftype

    return f


def convert_genbank(recs, seqid=None):
    """
    Converts BioPython SeqRecords obtained from a GENBANK file to a JSON representation.
    """
    global COUNTER

    # Start the counter from the one.
    reset_counter()

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

    logger.info("genbank conversion completed")
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
    name = f"{next(SEQNAME)}"
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
    attrs = dict(locus_tag=[name], start=start, end=end, type=ftype, strand=strand, location=location, operator=oper,
                 name=uid, id=uid)
    item[const.FEATURES] = [
        attrs
    ]
    data.append(item)

    return data


def json_view(params):
    """
    Prints json output to
    """
    for param in params:

        # Stop when data was not found.
        if not param.json:
            utils.error(f"data not found: {param.acc}")

        # Produce the full file when no parameters are set.
        if param.unset():
            text = json.dumps(param.json, indent=4)
            print(text)
        else:
            # Selects individual features.
            for item in param.json:
                feats = item[const.FEATURES]
                feats = filter_features(feats, param=param)
                text = json.dumps(list(feats), indent=4)
                print(text)


def parse_file(fname, seqid=None):
    """
    Parses a recognized file into a JSON representation
    """

    # Handle both compressed and uncompressed formats.
    stream = gzip.open(fname, 'rt') if fname.endswith(".gz") else open(fname, 'rt')

    # Detect extentions
    name, ext = os.path.splitext(fname)
    ext = ext.lower()

    # Split extension one more time if it looks like a compressed file.
    if ext == ".gz":
        name, ext = os.path.splitext(name)
        ext = ext.lower()

    # Cascade over the known file formats.
    if ext in (".gb", ".gbk", ".genbank"):
        recs = SeqIO.parse(stream, format=const.GENBANK)
        data = convert_genbank(recs, seqid=seqid)
    elif ext in (".fa", ".fasta"):
        recs = SeqIO.parse(stream, format=const.FASTA)
        data = convert_fasta(recs, seqid=seqid)
    else:
        utils.error(f"file format not recognized: {fname}")

    return data


if __name__ == "__main__":
    import doctest

    doctest.testmod()
