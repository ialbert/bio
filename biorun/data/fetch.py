"""
Fetches data from Entrez.
"""
import sys
import plac
from Bio import Entrez

Entrez.email = 'A.N.Other@example.com'


def error(msg):
    print(f"ERRROR: {msg}", file=sys.stderr)
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
def run(acc, db='nuccore', format='gb', mode='text'):
    stream = efetch(acc=acc, db=db, format=format, mode=mode)
    for line in stream:
        print (line, end="")
