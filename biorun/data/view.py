"""
Fetches data from Entrez.
"""
import sys
import os

import plac
from Bio import Entrez, SeqIO
from biorun import utils
from biorun.models import Sequence
from . import fetch

# The default logging function.
logger = utils.logger

# This email needs to be tunable.
Entrez.email = 'bio@bio.com'

def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def parse_genbank(stream):
    recs = SeqIO.parse(stream, utils.GENBANK)
    recs = map(lambda rec: Sequence(rec), recs)
    return recs



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


def process(acc, name='', fasta=False, gff=False):
    """
    Performs the processing of a single accession number.
    """

    # Open the stream to the data
    stream = fetch.get(acc=acc)

    if fasta:
        print_fasta(stream, name=name)
    elif gff:
        print_gff(stream, name=name)
    else:
        print_genbank(stream)

    return


@plac.pos('accs', "accession numbers")
@plac.opt('name', "name of the feature")
@plac.flg('fasta', "generate fasta file")
@plac.flg('gff', "generate a gff file")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(name='', fasta=False, gff=False, verbose=False, *accs):

    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    # The accession numbers may stored in files as well.
    collect = fetch.accs_or_file(accs)

    # Process each accession number.
    for acc in collect:
        process(acc, name=name, fasta=fasta, gff=gff)


