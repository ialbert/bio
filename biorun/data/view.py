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


def print_fasta(stream, name='', start=0, end=None, typ=None):
    """
    Prints the origin of a BioPython SeqRecord.
    """
    recs = parse_genbank(stream)

    # Print the origin for each record
    if not name:
        for item in recs:
            print(item.rec.format("fasta"))


def print_gff(stream, name='', start=0, end=None, typ=None):
    """
    Prints the origin of a BioPython SeqRecord.
    """
    recs = parse_genbank(stream)

    # Print the origin for each record
    for rec in recs:

        # Subselect by coordinates.
        feats = rec.features(start=start, end=end, name=name, typ=typ)

        # Generate the gff output
        for feat in feats:
            #print (feat)
            values = feat.as_gff(anchor=rec.id)
            values = map(str, values)
            print ("\t".join(values))



def print_genbank(stream):
    for line in stream:
        print(line, end="")


def process(acc, name='', fasta=False, gff=False, start=0, end=None, typ=''):
    """
    Performs the processing of a single accession number.
    """

    # Open the stream to the data
    cache, stream = fetch.get(acc=acc)

    if fasta:
        print_fasta(stream, name=name, start=start, end=end, typ=type)
    elif gff:
        print_gff(stream, name=name, start=start, end=end, typ=typ)
    else:
        print_genbank(stream)

    return


@plac.pos('acc', "accession numbers")
@plac.opt('name', "name of the feature")
@plac.opt('type', "the type of the feature")
@plac.opt('start', "start coordinate ")
@plac.opt('end', "end coordinate")
@plac.flg('fasta', "generate fasta file")
@plac.flg('gff', "generate a gff file")
@plac.flg('list', "lists the cached data")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(name='',  type='', start=0, end=0,  fasta=False, gff=False, verbose=False, list=False, *acc):

    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    # Produce the file listing
    if list:
        utils.print_file_list()

    start = start
    end = end or None
    # Process each accession number.
    for acx in acc:
        process(acx, name=name, fasta=fasta, gff=gff, start=start, end=end, typ=type)


