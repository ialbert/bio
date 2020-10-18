"""
Fetches data from Entrez and stores it in the local cache.
"""
import os, time, sys, re, json
import plac
from Bio import Entrez
from biorun import utils
import gzip
from biorun import models
from pprint import pprint
import collections

try:
    from Bio import SeqIO
except ImportError as exc:
    print(f"*** Error: {exc}", file=sys.stderr)
    print(f"*** Please install biopython: conda install -y biopython==1.76", file=sys.stderr)
    sys.exit(-1)

# The default logging function.
logger = utils.logger

Entrez.email = 'bio@bio.com'


def validate_ncbi(acc):
    """
   Catching a common error of passing a non existing file name as an accesion number.
    """
    base, ext = os.path.splitext(acc)
    if ext:
        try:
            # Version number must be integer
            int(ext[1:])
        except ValueError as exc:
            return False
    return True


@utils.time_it
def efetch(acc, db, format, mode='text'):
    if not validate_ncbi(acc):
        msg = f"invalid NCBI accession number: {acc}"
        utils.error(msg)

    try:
        logger.info(f"connecting to Entrez for {acc}")
        stream = Entrez.efetch(id=acc, db=db, rettype=format, retmode=mode)
        return stream
    except Exception as exc:
        msg = f"{exc} for efetch acc={acc} db={db} format={format} mode={mode}"
        utils.error(msg)

def gbk_to_json(gbk_name, json_name):
    """
    Transforms a GenBank file to a JSON file.
    """

    # Open GenBank file.
    inp_stream = gzip.open(gbk_name, 'rt') if gbk_name.endswith(".gz") else open(gbk_name, 'rt')

    # Convert genbank to a data structure.
    data = models.convert_genbank(inp_stream)

    # Save into a file.
    fp = gzip.open(json_name, 'wt', compresslevel=1)
    json.dump(data, fp)
    fp.close()

    logger.info(f"saved {json_name}")

def read_json_file(fname):
    stream = gzip.open(fname, 'rt')
    data = json.load(stream)
    stream.close()
    return data

def get_data(acc, db=None, format=utils.GENBANK, mode="text", update=False, rebuild=False):
    """
    Returns an open stream to the JSON file for a data.
    If the GenBank file does not exist it downloadsit as accession number from NCBI and converts
    it to json.
    """

    # Common mistake to pass an extra parameter.
    if acc.startswith("-"):
        msg = f"Invalid accession number: {acc}"
        utils.error(msg)

    # The accession number is a path to a GenBank file:
    if os.path.isfile(acc):
        logger.info(f"file {acc}")
        # Open the GenBank file.
        inp_stream = gzip.open(acc, 'rb') if acc.endswith(".gz") else open(acc, "rt")
        # Converts a GenBank to the internal JSON format the system uses.
        data = models.convert_genbank(inp_stream)
        return data

    # The JSON representation of the data.
    json_name = utils.resolve_fname(acc=acc , format="json")

    # GenBank representation of the data.
    gbk_name = utils.resolve_fname(acc=acc, format="gb")

    # Found the JSON representation of the file.
    if os.path.isfile(json_name) and not update and not rebuild:
        logger.info(f"found {json_name}")
        fp = gzip.open(json_name, 'rt')
        data = json.load(fp)
        return data

    # No JSON file but there is a genbank file. Convert, save the return the output.
    if os.path.isfile(gbk_name) and not update:
        logger.info(f"found {gbk_name}")
        gbk_to_json(gbk_name=gbk_name, json_name=json_name)
        data = read_json_file(json_name)
        return data

    # Accession numbers that are proteins.
    if acc[:3] in [ "AP_", "NP_", "YP_", "XP_", "WP_"]:
        db = "protein"
    else:
        db = db or "nuccore"

    # The data format.
    if format == utils.GENBANK:
        format = "gbwithparts"

    # Open the stream to the data.
    stream = efetch(acc=acc, db=db, format=format, mode=mode)

    # Save the strean to GenBank.
    utils.save_stream(stream=stream, fname=gbk_name)

    # Convert genbank to JSON.
    gbk_to_json(gbk_name=gbk_name, json_name=json_name)

    data = read_json_file(json_name)

    return data


def get_accessions(accs):
    """
    Returns a list that contains either the one accession or if any of these accessions is a
    valid filenames the content of those files.
    """
    collect = []
    for acc in accs:

        # No whitespaces allowed.
        acc = acc.strip()

        # Common mistake to pass an extra parameter to command line.
        if acc.startswith("-"):
            msg = f"Invalid accession number: {acc}"
            utils.error(msg)

        # Parse the file and extract accession numbers from it.
        if os.path.isfile(acc):

            # Get the lines (for sanity check limit at 100)
            lines = open(acc).readlines()[:100]

            # Remove whitespace from lines
            lines = [x.strip() for x in lines]

            # Remove comments from the file.
            lines = filter(lambda x: not x.startswith('#'), lines)

            # Remove empty lines from the file.
            lines = filter(None, lines)

            # Keep the first colums
            nums = [x.split()[0].strip() for x in lines]

            collect.extend(nums)
        else:
            collect.append(acc)

    return collect



@plac.pos('acc', "accession numbers")
@plac.opt('db', "database type", choices=["nuccore", "prot"])
@plac.flg('update', "download data again if it exists")
@plac.flg('quiet', "quiet mode, no output printed")
@plac.flg('build', "rebuilds the JSON representation")
def run(db='', alias='', update=False, quiet=False, build=False, *acc):

    # Set the verbosity level.
    utils.set_verbosity(logger, level=int(not quiet))

    # The accession numbers may stored in files as well.
    collect = get_accessions(acc)

    # Obtain the data for each accession number
    for acc in collect:

        data = get_data(acc=acc, db=db, update=update, rebuild=build)

        # A throttle to avoid accessing NCBI too quickly.
        if len(collect) > 1:
            time.sleep(1)
