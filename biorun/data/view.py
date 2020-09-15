"""
Fetches data from Entrez.
"""
import sys
import os

import plac
from Bio import Entrez
from biorun import utils

# The default logging function.
logger = utils.logger

Entrez.email = 'bio@example.com'

def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

@utils.time_it
def efetch(acc, db, format, mode='text'):
    try:
        stream = Entrez.efetch(id=acc, db=db, rettype=format, retmode=mode)
        return stream
    except Exception as exc:
        msg = f"{exc} for efetch acc={acc} db={db} format={format} mode={mode}"
        error(msg)

def get(acc, db='nuccore', format='gb', mode="text", update=False, stdout=False):
    """
    Downloads an accession number from NCBI to a file.
    Returns an open stream to the file.
    """
    # Resolve file based on the requested format.
    fname = utils.resolve_fname(acc=acc, format=format)

    # Override genbank mode to return the full file not just abbreviated.
    if format == 'gb':
        format = 'gbwithparts'

    # Return the file if the file already exists and update is not required.
    if os.path.isfile(fname) and not update:
        msg = f"file found at: {fname}"
        logger.info(msg)
        if stdout:
            for line in open(fname):
                print (line, end='')
    else:
        # Perform efetch and save the stream into the file.
        stream = efetch(acc=acc, db=db, format=format, mode=mode)
        utils.save_stream(stream=stream, fname=fname, stdout=stdout)

    stream = open(fname, 'rt')
    return stream

@plac.pos('acc', "accession number")
@plac.opt('db', "target database")
@plac.opt('format', "data format", choices=["gb", "fasta", "xml"])
@plac.flg('update', "update cached data if exists")
@plac.flg('verbose', "verbose mode, progress messages printed")
@plac.flg('prefetch', "saves data into cache, no output")
def run(acc, db='nuccore', format='gb', update=False, verbose=False, prefetch=False):

    utils.set_verbosity(logger, level=int(verbose))

    # Should it print output to the standard output
    stdout = not prefetch

    mode = "text"
    stream = get(acc=acc, db=db, format=format, mode=mode, update=update, stdout=stdout)

    stream.close()

    return
