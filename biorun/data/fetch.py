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


@plac.pos('acc', "accession number")
@plac.opt('db', "target database ")
@plac.opt('format', "output format")
@plac.opt('mode', "output mode")
@plac.opt('output_dir', "output directory")
@plac.flg('update', "overwrite existing data")
def run(acc, db='nuccore', format='gb', mode='text', output_dir=None, update=False):

    stream = efetch(acc=acc, db=db, format=format, mode=mode)

    # Resolve file name from accession number and download.
    outname = utils.resolve_fname(acc=acc, directory=output_dir)

    utils.download(stream=stream, outname=outname, overwrite=update)
