"""
Fetches data from Entrez.
"""
import sys
import os
import logging

import plac
from Bio import Entrez
from biorun import utils

Entrez.email = 'A.N.Other@example.com'

logger = logging.getLogger('bio')


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


def save_or_get(acc, db='nuccore', format='gb', mode='text', output_dir=None, output_name=None, update=False):
    # Resolve file name from accession number.
    outname = utils.resolve_fname(acc=acc, directory=output_dir, output_name=output_name)

    # Bail out if the output exists and update is False
    if os.path.isfile(outname) and not update:
        logger.info(f"File already exists at: {outname}")
        return outname

    # Preform efetch and return a generator object
    stream = efetch(acc=acc, db=db, format=format, mode=mode)

    # Download into outname.
    utils.save_file(stream=stream, outname=outname)

    return outname


@plac.pos('acc', "comma seperated list of accession numbers.")
@plac.opt('db', "target database ")
@plac.opt('format', "output format")
@plac.opt('mode', "output mode")
@plac.opt('output', "output directory")
@plac.opt('name', "output file name ( extension included )")
@plac.flg('overwrite', "overwrite existing data", abbrev='update')
def run(acc, db='nuccore', format='gb', mode='text', output=None, name=None, overwrite=False):

    output = save_or_get(acc=acc, db=db, format=format,
                         mode=mode, output_name=name,
                         output_dir=output,
                         update=overwrite)

    return
