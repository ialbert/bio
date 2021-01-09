"""
Handles functionality related to data storege.
"""
import sys, os, glob, re, gzip, json
from biorun import const, utils, objects, ncbi
from biorun.models import jsonrec
import biorun.libs.placlib as plac

# Module level logger.
logger = utils.logger

# A nicer error message on incorrect installation.
try:
    from Bio import SeqIO
except ImportError as exc:
    print(f"*** Error: {exc}", file=sys.stderr)
    print(f"*** This program requires biopython", file=sys.stderr)
    print(f"*** Install: conda install -y biopython>=1.78", file=sys.stderr)
    sys.exit(-1)

def resolve_fname(name, format='json'):
    """
    Resolve a file name given an accession number.
    """
    ext = format.lower()
    fname = f"{name}.{ext}.gz"
    fname = os.path.join(utils.DATADIR, fname)
    return fname


def delete_data(text):
    """
    Deletes data under a filename.
    """
    for name in text.split(","):
        fname = resolve_fname(name)
        if os.path.isfile(fname):
            os.remove(fname)
            logger.info(f"removed: {fname}")
        else:
            logger.info(f"file does not exist: {fname}")


def read_json_file(fname):
    """
    Returns the content of a JSON file.
    """
    fp = utils.gz_read(fname)
    data = json.load(fp)
    fp.close()
    return data


def save_json_file(fname, data):
    """
    Returns the content of a JSON file.
    """
    fp = utils.gz_write(fname)
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
        fp = utils.gz_write(json_name)
        json.dump(data, fp)
        fp.close()



def fetch_data(data, param):
    """
    Obtains data from NCBI. Fills each parameter with a json field.
    """

    db = "protein" if param.protein else "nuccore"

    for name in data:

        # Pretend no data if it is an update.
        json = None if param.update else get_json(name)

        # The data exists, nothing needs to be done.
        if json:
            continue

        # The JSON representation of the data.
        json_name = resolve_fname(name=name, format="json")

        # GenBank representation of the data.
        gbk_name = resolve_fname(name=name, format="gb")

        # Genome assembly data.
        if name.startswith("GCA") or name.startswith("GCF"):
            ncbi.genome(name=name, fname=gbk_name, update=param.update)
        else:
            # Genbank data.
            ncbi.genbank_save(name, db=db, fname=gbk_name)

        # Convert Genbank to JSON.
        data = jsonrec.parse_file(fname=gbk_name, seqid=param.seqid)

        # Save JSON file.
        save_json_file(fname=json_name, data=data)


def genbank_view(params):
    for param in params:
        altname = resolve_fname(param.acc, format="gb")

        if os.path.isfile(param.acc):
            stream = utils.gz_read(param.acc)
        elif os.path.isfile(altname):
            stream = utils.gz_read(altname)
        else:
            stream = []
            utils.error(f"data not found: {param.acc}")

        for line in stream:
            print(line, end='')


def get_json(name, seqid=None, inter=False, strict=False):
    """
    Attempts to return a JSON formatted data based on a name.
    """

    # Data is an existing path to a JSON file.
    if os.path.isfile(name):
        try:
            data = jsonrec.parse_file(name, seqid=seqid)
        except Exception as exc:
            logger.error(f"JSON parsing error for file {name}: {exc}")
            sys.exit(-1)
        return data

    # The JSON representation of the data.
    json_name = resolve_fname(name=name, format="json")

    # GenBank representation of the data.
    gbk_name = resolve_fname(name=name, format="gb")

    # Found the JSON representation of the file.
    if os.path.isfile(json_name):
        logger.info(f"found {json_name}")
        data = read_json_file(json_name)
        return data

    # There is no JSON file but there is a GenBank file.
    if os.path.isfile(gbk_name):
        logger.info(f"found {gbk_name}")
        data = jsonrec.parse_file(fname=gbk_name, seqid=seqid)
        data = save_json_file(fname=json_name, data=data)
        return data

    # Interactive input, make JSON from name
    if inter:
        data = jsonrec.make_jsonrec(name, seqid=seqid)
        return data

    # Raise error if in strict mode
    if strict:
        utils.error(f"data not found: {name}")
    return None


def rename_data(data, param, newname=None):
    """
    Rename data.
    """
    # It only makes sense to rename one data.
    newnames = newname.split(",")

    for name1, name2 in zip(data, newnames):
        src = resolve_fname(name=name1, format="json")
        dest = resolve_fname(name=name2, format="json")
        if os.path.isfile(src):
            logger.info(f"renamed {name1} as {name2}")
            os.rename(src, dest)
            if param.seqid:
                change_seqid(dest, seqid=param.seqid)
        else:
            logger.info(f"file not found: {src}")

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
        fsize = utils.human_size(os.path.getsize(path))
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


@plac.pos("data", "data names")
@plac.flg('fetch', "download data as accessions")
@plac.flg('update', "updates data in storage")
@plac.opt('rename', "rename the data")
@plac.opt('seqid', "set the sequence id of the data")
@plac.flg('protein', "use the protein database")
@plac.flg('verbose', "verbose mode")
def run(update=False, rename='', seqid='', protein=False, verbose=False, *data):
    """
    Fetches and manages data in storage.
    """


    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Reset counter (needed for consistency during testing).
    jsonrec.reset_counter()

    # A simple wrapper class to represent input parameters.
    param = objects.Param(seqid=seqid, rename=rename, start=1, protein=protein, update=update)

    # Fetch the data.
    fetch_data(data, param=param)

    # Renaming after fetching.
    if rename:
        rename_data(data, param=param, newname=rename)


@plac.opt('delete', "deletes foo from storage", metavar='foo')
@plac.flg('verbose', "verbose mode")
def data(delete, verbose=False):
    """
    Shows the data in the storage.

    Usage:

        bio data                   : lists the data
        bio data --delete foo      : deletes data called foo
        bio data --delete foo,bar  : deletes multiple datasets
    """
    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Reset counter (needed for consistency during testing).
    jsonrec.reset_counter()

    # Delete should be the first to execute.
    if delete:
        delete_data(delete)
    else:
        # Prints the data listing.
        print_data_list()