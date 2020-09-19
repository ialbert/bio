"""
Fetches data from Entrez and stores it in the local cache.
"""
import os, time, sys, re
import plac
from Bio import Entrez
from biorun import utils
import gzip

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


def get(acc, db='nuc', format=utils.GENBANK, mode="text", update=False, stdout=False):
    """
    Downloads an accession number from NCBI to a file.
    Returns an open stream to the file.
    """

    cache = True

    # Check if the accession number is a local file:
    if os.path.isfile(acc):
        logger.info(f"file {acc}")
        stream = gzip.open(acc, 'rb') if acc.endswith(".gz") else open(acc, "rt")
        return cache, stream

    if db == "nuc":
        db = "nuccore"
    else:
        db = "protein"

    if format == utils.GENBANK:
        format = "gbwithparts"

    # Flag indicating the source

    # Resolve file based on the requested format.
    fname = utils.resolve_fname(acc=acc, format=format)

    # Return the file if the file already exists and update is not required.
    if os.path.isfile(fname) and not update:
        logger.info(f"found {fname}")
    else:
        # Perform efetch and save the stream into the file.
        cache = False
        stream = efetch(acc=acc, db=db, format=format, mode=mode)
        utils.save_stream(stream=stream, fname=fname, stdout=stdout)

    stream = gzip.open(fname, 'rt')

    return cache, stream


def accs_or_file(accs):
    """
    Returns a list that contains either the accessions or if these accessions are valid filenames
    the content of those files.
    """
    collect = []
    for acc in accs:
        if os.path.isfile(acc):

            # Get the lines
            lines = open(acc).readlines()

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
            collect.append(acc.strip())

    return collect


@plac.pos('acc', "accession numbers")
@plac.opt('db', "database type", choices=["nuc", "prot"])
@plac.flg('update', "update cached data if exists")
@plac.flg('list', "lists the cached data")
@plac.flg('quiet', "quiet mode, no output printed")
def run(db='nuc', update=False, list=False, quiet=False, *acc):
    # Set the verbosity level.
    utils.set_verbosity(logger, level=int(not quiet))

    # Produce the file listing
    if list:
        utils.print_file_list()

    # The accession numbers may stored in files as well.
    collect = accs_or_file(acc)

    # Obtain the data for each accession number
    for acc in collect:

        cache, stream = get(acc=acc, db=db, format=utils.GENBANK, mode="text", update=update)

        # No need to keep the stream open.
        stream.close()

        # A throttle to avoid accessing NCBI too quickly.
        if not cache:
            time.sleep(1)
