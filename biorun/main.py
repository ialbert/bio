"""
The main job runner. Register functions here.
"""
import os, time, re
import sys
import plac
from biorun import VERSION
from biorun.const import *
from biorun import utils
from biorun.data import listing, view, storage

# Module level logger
logger = utils.logger


def smartname(text):
    """
    Splits an accession number by colon into acc:name
    """
    pass

@plac.flg('fasta', "produce FASTA format")
@plac.flg('gff', "produce GFF format", abbrev='G')
@plac.flg('fetch', "download data as accessions", abbrev='F')
@plac.flg('list', "list data in storage", abbrev='L')
@plac.flg('delete', "delete data in storage", abbrev='D')
@plac.flg('protein', "operate on proteins", abbrev='P')
@plac.flg('translate', "translate DNA to protein", abbrev='T')
@plac.opt('rename', "set the name", abbrev='R')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type")
@plac.opt('start', "start coordinate")
@plac.opt('end', "end coordinate")
@plac.opt('gene', "select features associated with gene" )
@plac.opt('match', "select features by rexep match")
@plac.flg('verbose', "verbose mode")
def base_runner(fasta=False, gff=False, fetch=False, protein=False, translate=False,
                delete=False, list=False, rename='',
                seqid='', start='', end='', type='', gene='', match='', verbose=False, *names):
    """
    bio - making bioinformatics fun again

    command line utility for manipulating bioinformatics data
    """

    # Check the names.
    names = storage.check_names(names)

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Delete the files from storage.
    if delete:
        storage.delete(names)

    # Get the data from Entrez.
    if fetch:
        db = "protein" if protein else "nuccore"
        storage.fetch(names, seqid=seqid, db=db)

    # The logging level may change within efetch to show progress.
    utils.set_verbosity(logger, level=int(verbose))

    # Renaming step before listing.
    if rename:
        storage.rename(names, seqid=seqid, newname=rename)

    # List the available data.
    if list:
        listing.print_data_list()


    # Populate the parameter list.
    param = utils.Param(start=start, end=end, seqid=seqid, protein=protein,
                        gff=gff, translate=translate, fasta=fasta, type=type, gene=gene, regexp=match)

    # Convert if no other command was given.
    convert = not(list or rename or delete or fetch)

    if convert:
        # Perform the data conversion
        view.convert_all(names, param=param)


def router():
    """
    Routes the tasks based on incoming parameters.
    """

    if ALIGN in sys.argv:

        # Delayed import to avoid missing library warning for other tasks.
        from biorun.align import pairwise

        # Drop the alignment request
        sys.argv.remove(ALIGN)

        # Add the help flag if otherwise empty.
        sys.argv += ["-h"] if len(sys.argv) == 1 else []

        # Call the pairwise aligner.
        plac.call(pairwise.run)

    else:

        # Add the help flag if otherwise empty.
        sys.argv += ["-h"] if len(sys.argv) == 1 else []
        plac.call(base_runner)


if __name__ == '__main__':
    router()
