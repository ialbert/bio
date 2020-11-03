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
        if os.path.isfile(fname):
            logger.info(f"removing: {fname}")
            os.remove(fname)
        else:
            logger.info(f"file does not exist: {fname}")


def check_names(names):
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
    fp = gzip.open(fname, 'rt') if fname.endswith(".gz") else open(fname, 'rt')
    data = json.load(fp)
    fp.close()
    return data

def save_json_file(fname, data):
    """
    Returns the content of a JSON file.
    """
    fp = gzip.open(fname, 'wt') if fname.endswith(".gz") else open(fname, 'wt')
    json.dump(data, fp)
    fp.close()
    logger.info(f"saved {fname}")
    return data


def change_seqid(json_name, seqid):
    """
    Changes the sequence id stored in a json file.
    """
    if os.path.isfile(json_name):
        data = read_json_file(json_name)
        for item in data:
            item[SEQID] = seqid
        fp = gzip.open(json_name, 'wt', compresslevel=1)
        json.dump(data, fp)
        fp.close()

def ncbi_efetch(name, gbk_name, db=None):
    """
    Connects to Entrez Direct to download data.
    """
    # Get the entire GenBank file.
    format, retmode = "gbwithparts", "text"

    # Guess accession numbers that are proteins.
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

    # Save the stream to GenBank.
    utils.save_stream(stream=stream, fname=gbk_name)


def fetch(names, seqid=None, db=None):
    """
    Obtains data from NCBI
    """
    # It is useful to see how large files are downloaded.
    utils.set_verbosity(logger, level=1)

    # Find names that do not exist
    names = filter(lambda n: not get_json(n), names)

    for name in names:
        # The JSON representation of the data.
        json_name = resolve_fname(name=name, format="json")

        # GenBank representation of the data.
        gbk_name = resolve_fname(name=name, format="gb")

        # Fetch and store genbank from remote site.
        ncbi_efetch(name, db=db, gbk_name=gbk_name)

        # Convert genbank to JSON.
        data = models.parse_file(fname=gbk_name, seqid=seqid)

        # Save JSON file.
        save_json_file(fname=json_name, data=data)

def get_json(name, seqid=None):
    """
    Attempts to return a JSON formatted data based on a name.
    """

    # Data is an existing path to a file.
    if os.path.isfile(name):
        data = models.parse_file(name, seqid=seqid)
        return data

    # Not a local file, attempt to resolve to storage.

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
        data = models.parse_file(fname=gbk_name, seqid=seqid)
        data = save_json_file(fname=json_name, data=data)
        return data

    return None

def rename(names, seqid=None, newname=None):

    # Empty list
    if not names:
        return

    # It only makes sense to rename one in case there are more.
    name = names[0]

    # We can only rename files that we have
    if get_json(name):
        src = resolve_fname(name=name, format="json")
        dest = resolve_fname(name=newname, format="json")
        if os.path.isfile(src):
            logger.info(f"moved {dest}")
            os.rename(src, dest)
            if seqid:
                change_seqid(dest, seqid=seqid)

        else:
            logger.error(f"not in storage: {src}")
    else:
        logger.error(f"not found: {name}")

