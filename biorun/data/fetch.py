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










def get_data(acc, db=None, format=utils.GENBANK, mode="text", update=False, rebuild=False, name=None):
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
    json_name = utils.resolve_fname(acc=(name or acc), format="json")

    # GenBank representation of the data.
    gbk_name = utils.resolve_fname(acc=acc, format="gb")

    # Found the JSON representation of the file.
    if os.path.isfile(json_name) and not update and not rebuild:
        logger.info(f"found {json_name}")
        data = read_json_file(json_name)
        return data

    # No JSON file but there is a genbank file. Convert, save the return the output.
    if os.path.isfile(gbk_name) and not update:
        logger.info(f"found {gbk_name}")
        gbk_to_json(gbk_name=gbk_name, json_name=json_name, seqid=name)
        data = read_json_file(json_name)
        return data



    # The data format.
    if format == utils.GENBANK:
        format = "gbwithparts"

    # Open the stream to the data.
    stream = efetch(acc=acc, db=db, format=format, mode=mode)

    # Save the strean to GenBank.
    utils.save_stream(stream=stream, fname=gbk_name)

    # Convert genbank to JSON.
    gbk_to_json(gbk_name=gbk_name, json_name=json_name, seqid=name)

    # Return the data
    data = read_json_file(json_name)

    return data


cs

@plac.pos('acc', "accession numbers")
@plac.opt('db', "database type", choices=["nuccore", "protein"])
@plac.opt('name', "rename the data (both the file and sequence id")
@plac.flg('update', "download data again if it exists")
@plac.flg('quiet', "quiet mode, no output printed")
@plac.flg('build', "rebuilds the JSON representation")
def run(db='', update=False, name='', quiet=False, build=False, *acc):
    # Set the verbosity level.
    utils.set_verbosity(logger, level=int(not quiet))

    # The accession numbers may stored in files as well.
    collect = get_accessions(acc)

    # Obtain the data for each accession number
    for acc in collect:

        data = get_data(acc=acc, db=db, update=update, rebuild=build, name=name)

        # A throttle to avoid accessing NCBI too quickly.
        if len(collect) > 1:
            time.sleep(1)
