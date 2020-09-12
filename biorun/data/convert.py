"""
Convert data from Entrez into different formats
"""
import os
import plac
from Bio import SeqIO
from biorun import utils, models
from biorun.data import fetch

logger = utils.logger

def converter(in_file, in_format, out_file, out_format):
    SeqIO.convert(in_file=in_file, in_format=in_format,
                  out_file=out_file,
                  out_format=out_format)
    return


def fasta_converter(in_file, in_format, out_file, verb=0):
    """
    Convert a given file with a format to
    """
    # Load and parse the file.
    input_handle = open(in_file, 'r')
    stream = SeqIO.parse(input_handle, in_format)

    # Make sure the file is in fasta format.
    formatter = lambda seq_record: f">{seq_record.id} {seq_record.description}\n{seq_record.seq}\n"

    # Download fasta file into out file
    utils.save_stream(stream, outname=out_file, verb=verb, formatter=formatter)

    return out_file


def print_origin(stream):
    """
    Prints the origin for each record in a BioPython seqrecord.
    """
    recs = models.parse_genbank(stream)
    for item in recs:
        print(item.rec.format("fasta"))

def convert(stream, outformat):

    pass

@plac.pos('acc', "accession number")
@plac.opt('db', "input target database ")
@plac.opt('format', "output data format")
@plac.flg('update', "update cached data if exists")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(acc, db='nuccore', format='fasta', verbose=False, update=False):

    utils.set_verbosity(logger, level=int(verbose))

    # Get a local copy of the genbank file or save_file it using Entrez.
    stream = fetch.get(acc=acc, db=db, format="gb",  update=update)

    print_origin(stream)
