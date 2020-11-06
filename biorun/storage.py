"""
Deals with the data storage.
"""
import sys, os, glob, re, gzip, json
from biorun import const, utils
from biorun.models import jsonrec

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


def delete(params):
    """
    Deletes data under a filename.
    """
    for p in params:
        fname = resolve_fname(p.name)
        if os.path.isfile(fname):
            logger.info(f"removing: {fname}")
            os.remove(fname)
        else:
            logger.info(f"file does not exist: {fname}")


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
            item[const.SEQID] = seqid
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


def fetch(params, seqid=None, db=None, update=False):
    """
    Obtains data from NCBI
    """

    for p in params:

        # Skip if exists (or not update).
        if p.json and not update:
            continue

        # The JSON representation of the data.
        json_name = resolve_fname(name=p.name, format="json")

        # GenBank representation of the data.
        gbk_name = resolve_fname(name=p.name, format="gb")

        # Fetch and store genbank from remote site.
        ncbi_efetch(p.name, db=db, gbk_name=gbk_name)

        # Convert genbank to JSON.
        data = jsonrec.parse_file(fname=gbk_name, seqid=seqid)

        # Save JSON file.
        save_json_file(fname=json_name, data=data)

def get_json(name, seqid=None, update=False, inter=False, strict=False):
    """
    Attempts to return a JSON formatted data based on a name.
    """

    # Data is an existing path to a file.
    if os.path.isfile(name):
        data = jsonrec.parse_file(name, seqid=seqid)
        return data

    # Not a local file, attempt to resolve to storage.

    # Report as not found if update is requested.
    if update:
        return None

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
        data = jsonrec.parse_file(fname=gbk_name, seqid=seqid)
        data = save_json_file(fname=json_name, data=data)
        return data

    # If not found and interactive mode create a JSON from the name itself.
    if inter:
        data = jsonrec.make_json(seq=name, seqid=seqid)
        return data

    # At this point the data was not found
    if strict:
        utils.error(f"data not found: {name}")

    return None

def rename(params, seqid=None, newname=None):

    # Empty list
    if not params:
        return

    # It only makes sense to rename one of the many
    name = params[0].name

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


def print_data_list():
    """
    Returns a list of the files in the data directory
    """
    pattern = os.path.join(os.path.join(utils.DATADIR, '*.json.gz'))
    matched = glob.glob(pattern)

    # Extract the definition from the JSON without parsing it.
    patt = re.compile(r'(definition\":\s*)(?P<value>\".+?\")')
    collect = []

    for path in matched:
        fsize = utils.sizeof_fmt(os.path.getsize(path))
        base, fname = os.path.split(path)
        fname = fname.rsplit(".", maxsplit=2)[0]

        # Parse the first N lines
        stream = gzip.open(path, 'rt') if path.endswith('gz') else open(path, 'rt')
        text = stream.read(1000)
        match = patt.search(text)

        title = match.group("value") if match else ''
        title = title.strip('", ')

        # Trim the title
        stitle = title[:100]
        stitle = stitle + "..." if len(title) != len(stitle) else stitle

        collect.append((str(fsize), f"{fname:10s}", stitle))

    collect = sorted(collect, key=lambda x: x[2])
    for row in collect:
        line = "\t".join(row)
        print(line)
