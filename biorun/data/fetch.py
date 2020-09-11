"""
Fetches data from Entrez.
"""
import sys
import os

import plac
from Bio import Entrez
from biorun import utils

Entrez.email = 'A.N.Other@example.com'


def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def efetch(acc, db, format, mode='text'):
    try:

        stream = Entrez.efetch(id=acc, db=db, rettype=format, retmode=mode)
        for line in stream:
            yield line

    except Exception as exc:
        msg = f"{exc} for efetch acc={acc} db={db} format={format} mode={mode}"
        error(msg)


def save_or_get(acc, db='nuccore', format='gb', mode='text', output_dir=None, update=False, verbosity=0):
    # Resolve file name from accession number.
    outname = utils.resolve_fname(acc=acc, directory=output_dir)

    # Bail out if the output exists and update is False
    if os.path.isfile(outname) and not update:
        msg = f"File already exists at: {outname}"
        utils.print_message(msg=msg, styles=[utils.BOLD], verb=verbosity)
        return outname

    # Preform efetch and return a generator object
    stream = efetch(acc=acc, db=db, format=format, mode=mode)

    # Download into outname.
    utils.save_file(stream=stream, outname=outname, verb=verbosity)

    return outname


@plac.pos('acc', "accession number")
@plac.opt('db', "target database ")
@plac.opt('format', "output format")
@plac.opt('mode', "output mode")
@plac.opt('output_dir', "output directory")
@plac.flg('update', "overwrite existing data")
@plac.opt('verbosity', "verbosity level", type=int)
def run(acc, db='nuccore', format='gb', mode='text', output_dir=None, update=False, verbosity=0):

    output = save_or_get(acc=acc, db=db, format=format,
                         mode=mode,
                         output_dir=output_dir,
                         update=update,
                         verbosity=verbosity)

    return
