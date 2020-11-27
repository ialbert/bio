"""
Deals with converting between formats.
"""
import biorun.libs.placlib as plac

from biorun import utils, const, storage, objects
from biorun.models import fastarec, gffrec, jsonrec

# Module level logger
logger = utils.logger


@plac.pos("acc", "data names")
@plac.flg('fasta', "produce FASTA format", abbrev='F')
@plac.flg('gff', "produce GFF format", abbrev='G')
@plac.flg('fetch', "download data as accessions")
@plac.flg('list', "list data in storage")
@plac.flg('delete', "delete data in storage")
@plac.flg('update', "updates data in storage")
@plac.flg('protein', "operate on proteins", abbrev='P')
@plac.flg('translate', "translate DNA to protein", abbrev='T')
@plac.flg('transcribe', "transrcribe DNA to RNA", abbrev='X')
@plac.opt('rename', "set the name", abbrev='r')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type")
@plac.opt('start', "start coordinate")
@plac.opt('end', "end coordinate")
@plac.opt('gene', "select features associated with gene")
@plac.opt('name', "select features associated with a name", abbrev='N')
@plac.opt('match', "select features by rexep match")
@plac.flg('inter', "interactive (data from command line)", abbrev='i')
@plac.flg('origin', "use the origin (source) sequences", abbrev='O')
@plac.flg('reverse', "reverse sequence", abbrev='E')
@plac.flg('genbank', "show the genbank file if exists", abbrev='K')
@plac.flg('complement', "complement sequence", abbrev='C')
@plac.flg('revcomp', "reverse complement sequence", abbrev='R')
@plac.flg('verbose', "verbose mode")
def run(fasta=False, gff=False, genbank=False, fetch=False, update=False, delete=False, list=False, protein=False,
              translate=False, transcribe=False,
              reverse=False, complement=False, revcomp=False, rename='', seqid='', start='', end='', type='', gene='', name='',
              match='', inter=False, origin=False,
              verbose=False, *acc):
    """
    bio - making bioinformatics fun again

    Command line utilities for manipulating bioinformatics data.

    Subcommands with additional help:

        bio --align -h
        bio --taxon -h
        bio --sra -h

    """

    def make_param(acc):
        """
        Creates a parameter for each accession.

        """
        # Set the verbosity
        utils.set_verbosity(logger, level=int(verbose))

        # A very common error to pass a fragment as
        if acc.startswith("-"):
            msg = f"Invalid accession number: {acc}"
            utils.error(msg)

        # A simple wrapper class to carry all parameters around.
        p = objects.Param(start=start, end=end, seqid=seqid, protein=protein, revcomp=revcomp,
                          update=update, acc=acc, gff=gff, translate=translate, reverse=reverse,
                          complement=complement, origin=origin, name=name,
                          fasta=fasta, type=type, gene=gene, regexp=match, transcribe=transcribe)

        # Fill the json data for the parameter.
        p.json = storage.get_json(p.acc, seqid=seqid, inter=inter)
        return p

    # Allow commas in numbers, or sizes like 10Kb
    start = utils.parse_number(start)
    end = utils.parse_number(end)

    # Make a list of parameters for each name.
    params = [make_param(a) for a in acc]

    # Delete should be the first to execute.
    if delete:
        storage.delete(params)

    # Fetch to be perfomed before renaming.
    if fetch:
        db = "protein" if protein else "nuccore"
        storage.fetch(params, seqid=seqid, db=db, update=update)

    # Renaming needs to be performed before listing.
    if rename:
        storage.rename(params, seqid=seqid, newname=rename)

    # List the available data.
    if list:
        storage.print_data_list()

    # Stop after storage related actions
    if (list or rename or delete):
        return

    # Sequence operation
    seqop = reverse or complement or revcomp

    # Decide which type of conversion based on incoming parameters.
    if genbank:
        storage.genbank_view(params)
    elif origin or fasta or protein or translate or seqop:
        fastarec.fasta_view(params)
    elif gff:
        gffrec.gff_view(params)
    elif not fetch:
        jsonrec.json_view(params)
