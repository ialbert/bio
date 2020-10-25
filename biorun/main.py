"""
The main job runner. Register functions here.
"""
import os, time
import sys
import plac
from biorun import VERSION
from biorun import utils
from biorun.data import listing, view, storage

# Module level logger
logger = utils.logger

def zero_based(start, end):
    """
    Shift to zero based coordinate system.
    """

    # Don't allow zero for start.
    if start == 0:
        utils.error(f"start={start} may not be  zero")

    # Default value for start
    start = int(start) if start else 1

    # Default value for end.
    end = None if not end else int(end)

    # Shift the coordinates
    start = start - 1 if start > 0 else start
    return start, end


class Param(object):
    """
    Parameter representation that have grown too numerous to pass individually.
    """

    def __init__(self, **kwds):
        self.start = self.end = self.seqid = None
        self.gff = self.protein = self.fasta = self.translate = None
        self.name = self.gene = self.type = self.regexp = None
        self.__dict__.update(kwds)
        self.start, self.end = zero_based(start=self.start, end=self.end)

    def unset(self):
        """
        Feature filtering parameters not set.
        """
        return not (self.start or self.end or self.type or self.gene or self.regexp)

    def __str__(self):
        return str(self.__dict__)

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
@plac.flg('store', "places a file into storage", abbrev='X')
@plac.opt('rename', "set the name", abbrev='R')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type")
@plac.opt('start', "start coordinate")
@plac.opt('end', "end coordinate")
@plac.opt('gene', "select features associated with gene" )
@plac.opt('match', "select features by rexep match")
@plac.opt('align', "alignment mode", choices=['global', 'local'])
@plac.flg('verbose', "verbose mode")
def run(fasta=False, gff=False, fetch=False, protein=False, translate=False,
        delete=False,  list=False, store=False, rename='',
        seqid='', start='', end='', type='', gene='', match='', align='', verbose=False, *names):
    """
    bio - making bioinformatics fun again

    command line utility for manipulating bioinformatics data
    """

    # Check the names.
    names = storage.validate_names(names)

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Delete the files from storage.
    if delete:
        storage.delete(names)

    # Get the data from Entrez.
    if fetch:
        db = "protein" if protein else None
        storage.fetch(names, seqid=seqid, db=db)

    # The logging version may change within efetch to show progress.
    utils.set_verbosity(logger, level=int(verbose))

    # Renaming step.
    if rename:
        storage.rename(names, seqid=seqid, newname=rename)

    # List the available data.
    if list:
        listing.print_data_list()
        sys.exit()

    # Stop here.
    if list or rename:
        sys.exit()

    # Populate the parameter list.
    param = Param(start=start, end=end, seqid=seqid, protein=protein,
                        gff=gff, translate=translate, fasta=fasta, type=type, gene=gene, regexp=match)

    # Decide if it is a data conversion
    convert = not(align or delete)

    if convert:
        # Perform the data conversion
        view.convert_all(names, param=param)


def toplevel():
    """
    Runs the toplevel function.
    """
    plac.call(run)


if __name__ == '__main__':
    toplevel()
