"""
The main job runner. Register functions here.
"""
import os, time, re
import sys
import plac

from biorun import utils, const
from biorun.data import listing, storage
from biorun.data import fastarec, gffrec, jsonrec

# Module level logger
logger = utils.logger


@plac.flg('fasta', "produce FASTA format")
@plac.flg('gff', "produce GFF format", abbrev='G')
@plac.flg('fetch', "download data as accessions", abbrev='F')
@plac.flg('list', "list data in storage", abbrev='L')
@plac.flg('delete', "delete data in storage", abbrev='D')
@plac.flg('update', "updates data in storage", abbrev='U')
@plac.flg('protein', "operate on proteins", abbrev='P')
@plac.flg('translate', "translate DNA to protein", abbrev='T')
@plac.opt('rename', "set the name", abbrev='R')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type")
@plac.opt('start', "start coordinate")
@plac.opt('end', "end coordinate")
@plac.opt('gene', "select features associated with gene")
@plac.opt('match', "select features by rexep match")
@plac.flg('verbose', "verbose mode")
def converter(fasta=False, gff=False, fetch=False, update=False, protein=False, translate=False,
              delete=False, list=False, rename='', seqid='', start='', end='', type='', gene='', match='',
              verbose=False, *acc):
    """
    bio - making bioinformatics fun again

    command line utility for manipulating bioinformatics data
    """

    def make_param(name):
        """
        Creates a parameter for each accession.

        """
        # A very common error. Catch it here.
        if name.startswith("-"):
            msg = f"Invalid accession number: {name}"
            utils.error(msg)

        # A simple wrapper class to carry all parameters around.
        p = utils.Param(start=start, end=end, seqid=seqid, protein=protein,
                        update=update, name=name, gff=gff, translate=translate,
                        fasta=fasta, type=type, gene=gene, regexp=match)

        # Fill the json data for the name.
        p.json = storage.get_json(name, seqid=seqid)
        return p

    # Make a list of parameters for each name.
    params = [make_param(n) for n in acc]

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Delete the files from storage.
    if delete:
        storage.delete(params)

    # Get the data from Entrez.
    if fetch:
        db = "protein" if protein else "nuccore"
        storage.fetch(params, seqid=seqid, db=db, update=update)

    # Renaming step before listing.
    if rename:
        storage.rename(params, seqid=seqid, newname=rename)

    # List the available data.
    if list:
        listing.print_data_list()

    # Condition to exit.
    exit = (list or rename or delete or fetch)

    if exit:
        return

    # Looks like  conversion
    if fasta:
        fastarec.fasta_view(params)
    elif gff:
        gffrec.gff_view(params)
    else:
        jsonrec.json_view(params)


def router():
    """
    Routes the tasks based on incoming parameters.
    """

    # Alignment requested.
    if const.ALIGN in sys.argv:

        # Drop the alignment command from paramters.
        sys.argv.remove(const.ALIGN)

        # Delayed import to avoid missing library warning for other tasks.
        from biorun.align import pairwise

        # Add the help flag if otherwise empty.
        sys.argv += ["-h"] if len(sys.argv) == 1 else []

        # Call the pairwise aligner.
        plac.call(pairwise.run)

    # Default action is to convert a file.
    else:

        # Add the help flag if otherwise empty.
        sys.argv += ["-h"] if len(sys.argv) == 1 else []
        plac.call(converter)


if __name__ == '__main__':
    router()
