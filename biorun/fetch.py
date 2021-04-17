"""
Handles functionality related to data storege.
"""
import sys, os
from biorun import utils
import biorun.libs.placlib as plac
from itertools import *

# Module level logger.
logger = utils.logger

# A nicer error message on incorrect installation.
try:
    from Bio import Entrez
except ImportError as exc:
    print(f"*** Error: {exc}", file=sys.stderr)
    print(f"*** This program requires biopython", file=sys.stderr)
    print(f"*** Install: conda install -y biopython>=1.78", file=sys.stderr)
    sys.exit(-1)

# This is to silence Biopython warning.
Entrez.email = 'not set'

def info(msg='', quiet=False, end="\n", file=sys.stderr):
    if not quiet:
        print(f"*** {msg}", file=file, end=end)


def ncbi_efetch(ids, db, rettype='gbwithparths', retmode='text', quiet=False):
    """
    Connects to Entrez Direct to download data.
    """

    # Guess accession numbers that are proteins.
    # https: // www.ncbi.nlm.nih.gov /Sequin /acc.html

    size = 2000

    idtext = ",".join(ids)

    try:
        info(f"connecting to Entrez for {len(ids)} records", quiet=quiet)
        stream = Entrez.epost(db=db, id=idtext)
        results = Entrez.read(stream)
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
        stream = Entrez.efetch(webenv=webenv, query_key=query_key, db=db, rettype=rettype, retmode=retmode)
    except KeyError as exc:
        msg1 = f"{exc}"
        msg2 = f"efetch acc={ids} db={db} rettype={rettype} retmode={retmode}"
        logger.error(msg1)
        logger.error(msg2)

        sys.exit()

    recnum = 0
    step = count(1)
    for index, line in zip(step, stream):

        # Print line to standard output.
        print(line, end="")

        # New record was hit.
        start = line.startswith("LOCUS")

        # Increment record number.
        if start:
            recnum += 1

        # Download progress messages.
        if start or (index % size) == 0:
            info(f"{recnum} records, processed {index:,d} lines\r", quiet=quiet, end='')

    # Final message
    info(f"downloaded {recnum} records, processed {index:,d} lines\r", quiet=quiet, end='')

    # Move cursor on new line.
    info(quiet=quiet)


# Guess accession numbers that are proteins.
# https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/

def validator(acc):
    '''
    Validate accession numbers
    '''
    # Get rid of version number
    acc = acc.split(".")[0]

    # Last 5 digits need to be numerals
    size = len(acc)

    try:
        int(acc[-5:])
    except ValueError as exc:
        return False

    return size >= 5


@plac.opt('db', "database name")
@plac.opt('chunk', "chunk size (max records per connection)")
@plac.opt('format_', "output format")
@plac.flg('quiet', "quiet mode")
def run(db='nuccore', format_="gb", quiet=False, chunk=1000, *data):
    """
    Fetches and manages data in storage.
    """

    ids = []
    for elem in data:
        if os.path.isfile(elem):
            elems = map(lambda x: x.strip().split()[0], open(elem))
            ids.extend(elems)
        else:
            ids.append(elem)

    ids = filter(validator, ids)
    ids = list(ids)

    chunks = [ids[x:x + chunk] for x in range(0, len(ids), chunk)]

    if len(chunks) > 1:
        info(f"input data contains {len(ids)} records", quiet=quiet)

    # A shortcut notation
    format_ = 'gbwithparts' if format_ == 'gb' else format_
    for subset in chunks:
        ncbi_efetch(subset, db=db, rettype=format_, quiet=quiet)
