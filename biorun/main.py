"""
The main job runner. Register functions here.
"""
import sys
import plac

from biorun import utils, const, storage, objects
from biorun.models import fastarec, gffrec, jsonrec

# Module level logger
logger = utils.logger


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
@plac.opt('match', "select features by rexep match")
@plac.flg('inter', "interactive (data from command line)", abbrev='i')
@plac.flg('reverse', "reverse sequence", abbrev='E')
@plac.flg('complement', "complement sequence", abbrev='C')
@plac.flg('revcomp', "reverse complement sequence", abbrev='R')
@plac.flg('verbose', "verbose mode")
def converter(fasta=False, gff=False, fetch=False, update=False, delete=False, list=False, protein=False,
              translate=False, transcribe=False,
              reverse=False, complement=False, revcomp=False, rename='',  seqid='', start='', end='', type='', gene='', match='', inter=False,
              verbose=False, *acc):
    """
    bio - making bioinformatics fun again

    command line utility for manipulating bioinformatics data
    """

    def make_param(name):
        """
        Creates a parameter for each accession.

        """
        # Set the verbosity
        utils.set_verbosity(logger, level=int(verbose))

        # A very common error to pass a fragment as
        if name.startswith("-"):
            msg = f"Invalid accession number: {name}"
            utils.error(msg)

        # A simple wrapper class to carry all parameters around.
        p = objects.Param(start=start, end=end, seqid=seqid, protein=protein, revcomp=revcomp,
                        update=update, name=name, gff=gff, translate=translate, reverse=reverse, complement=complement,
                        fasta=fasta, type=type, gene=gene, regexp=match, transcribe=transcribe)

        # Fill the json data for the parameter.
        p.json = storage.get_json(p.name, seqid=seqid, inter=inter)
        return p

    # Make a list of parameters for each name.
    params = [make_param(n) for n in acc]

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

    # Decide which type of conversion based on incoming parameters.
    if fasta or protein or translate:
        fastarec.fasta_view(params)
    elif gff:
        gffrec.gff_view(params)
    elif not fetch:
        jsonrec.json_view(params)

def proofreader(value):
    """
    Bridges the gap between short and longforms. Allows the use of both by remapping to canonical forms.
    Remaps -start to --start, --F to -F
    """

    # Longform parameter
    long = value.startswith("--")

    # Shortform parameter.
    short = not long and value.startswith("-")

    # Is it one character.
    onechar = len(value.strip("-")) == 1

    try:
        # Can the value be converted to a number
        float(value)
        isnum = True
    except Exception as exc:
        isnum = False

    # Single character but long form. Drop a leading dash.
    if onechar and long and not isnum:
        value = value[1:]

    # Short form but more than one char. Add a dash.
    if short and not onechar and not isnum:
        value = f"-{value}"

    return value

def router():
    """
    Routes the tasks based on incoming parameters.
    """

    # Allow multiple forms of parameters to be used.
    sys.argv = list(map(proofreader, sys.argv))

    # Alignment requested.
    if const.ALIGN in sys.argv:

        # Drop the alignment command from paramters.
        sys.argv.remove(const.ALIGN)

        # Delayed import to allow other functionality to work even when the parasail library is missing.
        from biorun.methods import align

        # Add the help flag if no other information is present.
        if len(sys.argv) == 1:
            sys.argv.append("-h")

        # Delegate parameter parsing to aligner.
        plac.call(align.run)

    elif  const.TAXON in sys.argv:
        from biorun.models import taxdb

        sys.argv.remove(const.TAXON)

        plac.call(taxdb.run)

    # Default action is to convert a file.
    else:

        # Add the help flag if no other information is present.
        if len(sys.argv) == 1:
            sys.argv.append("-h")

        # Delegate parameter parsing to converter.
        plac.call(converter)


if __name__ == '__main__':
    router()
