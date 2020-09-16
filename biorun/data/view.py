"""
Fetches data from Entrez.
"""
import sys
import os

import plac
from Bio import Entrez, SeqIO
from biorun import utils
from biorun.models import Sequence

# The default logging function.
logger = utils.logger

Entrez.email = 'bio@example.com'

def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def parse_genbank(stream):
    recs = SeqIO.parse(stream, utils.GENBANK)
    recs = map(lambda rec: Sequence(rec), recs)
    return recs


@utils.time_it
def efetch(acc, db, format, mode='text'):
    try:
        logger.info("connecting to edirect")
        stream = Entrez.efetch(id=acc, db=db, rettype=format, retmode=mode)
        logger.info("returning from edirect")
        return stream
    except Exception as exc:
        msg = f"{exc} for efetch acc={acc} db={db} format={format} mode={mode}"
        error(msg)


def print_fasta(stream, name=''):
    """
    Prints the origin of a BioPython SeqRecord.
    """
    recs = parse_genbank(stream)

    # Print the origin for each record
    if not name:
        for item in recs:
            print(item.rec.format("fasta"))

def print_gff(stream, name=''):
    """
    Prints the origin of a BioPython SeqRecord.
    """
    recs = parse_genbank(stream)

    # Print the origin for each record
    for rec in recs:
        for ival in rec:
            values = ival.data.as_gff(anchor=rec.name)
            values = map(str, values)
            print ("\t".join(values))

def print_genbank(stream):
    for line in stream:
        print(line, end="")


def get(acc, db='nuccore', format=utils.GENBANK, mode="text", update=False, stdout=False):
    """
    Downloads an accession number from NCBI to a file.
    Returns an open stream to the file.
    """

    if format == utils.GENBANK:
        format = "gbwithparts"

    # Resolve file based on the requested format.
    fname = utils.resolve_fname(acc=acc, format=format)


    # Return the file if the file already exists and update is not required.
    if os.path.isfile(fname) and not update:
        msg = f"file found at: {fname}"
        logger.info(msg)
    else:
        # Perform efetch and save the stream into the file.
        stream = efetch(acc=acc, db=db, format=format, mode=mode)
        utils.save_stream(stream=stream, fname=fname, stdout=stdout)

    stream = open(fname, 'rt')
    return stream

@plac.pos('acc', "accession number")
@plac.opt('db', "target database")
@plac.opt('name', "name of the feature")
@plac.flg('fasta', "generate fasta file")
@plac.flg('gff', "generate a gff file")
@plac.flg('update', "update cached data if exists")
@plac.flg('verbose', "verbose mode, progress messages printed")
@plac.flg('prefetch', "saves data into cache, no output")
def run(acc, db='nuccore', name='', fasta=False, gff=False,  update=False, verbose=False, prefetch=False):

    utils.set_verbosity(logger, level=int(verbose))

    parts = acc.split(":")
    if len(parts) == 2:
        acc = parts[0]
        name = name or parts[1]


    mode = "text"
    format = utils.GENBANK

    stream = get(acc=acc, db=db, format=format, mode=mode, update=update)

    # No further action taken.
    if prefetch:
        return

    if fasta :
        print_fasta(stream, name=name)
    elif gff:
        print_gff(stream, name=name)
    else:
        print_genbank(stream)

    return
