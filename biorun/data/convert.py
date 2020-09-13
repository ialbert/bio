"""
Convert data from Entrez into different formats
"""
import os
import plac
from Bio import SeqIO
from biorun import utils, models
from biorun.data import fetch

logger = utils.logger

def print_origin(stream):
    """
    Prints the origin for each record in a BioPython seqrecord.
    """
    recs = models.parse_genbank(stream)
    for item in recs:
        print(item.rec.format("fasta"))

def print_gff(stream):
    recs = models.parse_genbank(stream)
    for item in recs:
        for interval in item:
            print (interval)


@plac.pos('acc', "accession number")
@plac.opt('db', "input target database ")
@plac.opt('format', "output format", choices=[utils.FASTA, utils.GFF])
@plac.flg('update', "update cached data if exists")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(acc, db='nuccore', format='fasta', verbose=False, update=False):

    utils.set_verbosity(logger, level=int(verbose))

    # Get a local copy of the genbank file or save_file it using Entrez.
    stream = fetch.get(acc=acc, db=db, format="gb",  update=update)

    if format == utils.FASTA:
        print_origin(stream)
    elif format == utils.GFF:
        print_gff(stream)
