"""
The main job runner. Register functions here.
"""
import os, time
import sys
import plac
from biorun import VERSION
from biorun import utils
from biorun.data import listing
from biorun.data import storage

# Module level logger
logger = utils.logger

@plac.flg('fasta', "produce FASTA format")
@plac.flg('gff', "produce GFF format", abbrev='G')
@plac.flg('fetch', "download data as accessions", abbrev='F')
@plac.flg('list', "list data in storage", abbrev='L')
@plac.flg('delete', "delete data in storage", abbrev='D')
@plac.flg('protein', "operate on proteins", abbrev='P')
@plac.flg('translate', "translate DNA to protein", abbrev='T')
@plac.flg('store', "places a file into storage", abbrev='X')
@plac.opt('rename', "rename data in storage", abbrev='R')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type")
@plac.opt('start', "start coordinate")
@plac.opt('end', "end coordinate")
@plac.opt('gene', "select features associated with gene" )
@plac.opt('match', "select features by rexep match")
@plac.opt('align', "alignment mode", choices=['global', 'local'])
def run(fasta=False, gff=False, fetch=False, protein=False, translate=False,
        delete=False,  list=False, store=False, rename='',
        seqid='', start='', end='', type='', gene='', match='', align='', *names):
    """
    bio - making bioinformatics fun again

    command line utility for manipulating bioinformatics data
    """

    # Check the names.
    names = storage.validate_names(names)

    # Delete the files from storage.
    if delete:
        storage.delete(names)

    # Get the data from Entrez.
    if fetch:
        utils.set_verbosity(logger, level=1)
        db = "protein" if protein else None
        storage.fetch(names, seqid=seqid, db=db)

    if rename:
        storage.rename(names, seqid=seqid, newname=rename)

    # Listing has the higher priority.
    if list:
        listing.print_data_list()

    # Perform the


def toplevel():
    """
    Runs the toplevel function.
    """
    plac.call(run)


if __name__ == '__main__':
    toplevel()
