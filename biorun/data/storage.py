"""
Deals with the data storage.
"""
import sys, os, gzip, json
from biorun import utils, models
from biorun.const import *

# Module level logger.
logger = utils.logger

try:
    from Bio import SeqIO
    from Bio import Entrez
except ImportError as exc:
    print(f"*** Error: {exc}", file=sys.stderr)
    print(f"*** Please install biopython: conda install -y biopython==1.76", file=sys.stderr)
    sys.exit(-1)

Entrez.email = 'bio@bio.com'

def resolve_fname(name, format='json'):
    """
    Resolve a file name given an accession number.
    """
    ext = format.lower()
    fname = f"{name}.{ext}.gz"
    fname = os.path.join(utils.DATADIR, fname)
    return fname


def delete(names):
    """
    Deletes data under a filename.
    """
    for name in names:
        fname = resolve_fname(name)
        logger.info(f"removing: {fname}")
        if os.path.isfile(fname):
            os.remove(fname)


def validate_names(names):
    """
    Catches common user errors.
    """
    for name in names:
        if name.startswith("-"):
            msg = f"Invalid accession number: {name}"
            utils.error(msg)

    return names


def read_json_file(fname):
    """
    Returns the content of a JSON file.
    """
    stream = gzip.open(fname, 'rt') if fname.endswith(".gz") else open(fname, 'rt')
    data = json.load(stream)
    stream.close()
    return data


def gbk_to_json(gbk_name, json_name, seqid=None):
    """
    Transforms a GenBank file to a JSON file.
    """

    # Open GenBank file.
    inp_stream = gzip.open(gbk_name, 'rt') if gbk_name.endswith(".gz") else open(gbk_name, 'rt')

    # Convert genbank to a data structure.
    data = models.convert_genbank(inp_stream, seqid=seqid)

    # Save into a file.
    fp = gzip.open(json_name, 'wt', compresslevel=1)
    json.dump(data, fp)
    fp.close()

    logger.info(f"saved {json_name}")

def efetch(name, db=None):
    """
    Connects to Entrez Direct to download data.
    """
    # Get the entire GenBank file.
    format, retmode = "gbwithparts", "text"

    # Accession numbers that are proteins.
    if name[:2] in ["AP", "NP", "YP", "XP", "WP", "AK"]:
        db = db or "protein"
    else:
        db = db or "nuccore"

    try:
        logger.info(f"connecting to Entrez for {name}")
        stream = Entrez.efetch(id=name, db=db, rettype=format, retmode=retmode)
    except Exception as exc:
        msg = f"{exc} for efetch acc={name} db={db} format={format} mode={retmode}"
        utils.error(msg)

    # The name of the file to save it locally to.
    gbk_name = resolve_fname(name=name, format="gb")

    # Save the stream to GenBank.
    utils.save_stream(stream=stream, fname=gbk_name)

    return gbk_name

def fetch(names, seqid=None, db=None):

    # Find names that do not exist
    names = filter(lambda n: not get_json(n), names)

    for name in names:
        # The JSON representation of the data.
        json_name = resolve_fname(name=name, format="json")

        # Need to fetch genbank from remote site.
        gbk_name = efetch(name, db=db)

        # Copy GenBank to JSON
        gbk_to_json(gbk_name=gbk_name, json_name=json_name, seqid=seqid)


def get_json(name, seqid=None):
    """
    Attempts to return a JSON formatted data for a name.
    """

    # The JSON representation of the data.
    json_name = resolve_fname(name=name, format="json")

    # GenBank representation of the data.
    gbk_name = resolve_fname(name=name, format="gb")

    # Found the JSON representation of the file.
    if os.path.isfile(json_name):
        logger.info(f"found {json_name}")
        data = read_json_file(json_name)
        return data

    # No JSON file but there is a genbank file.
    if os.path.isfile(gbk_name):
        logger.info(f"found {gbk_name}")
        gbk_to_json(gbk_name=gbk_name, json_name=json_name, seqid=seqid)
        data = read_json_file(json_name)
        return data

    return None

def rename(names, seqid=None, newname=None):

    # Renames only the first
    pass

