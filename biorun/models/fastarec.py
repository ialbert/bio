"""
Handles FASTA related outputs
"""

import biorun.libs.placlib as plac
from biorun.models import jsonrec
from biorun import utils, const, storage, objects

logger = utils.logger


def get_fasta(item, param):
    """
    Returns a fasta seqrecord for a JSON item.
    """

    # Ignore translation warnings when the sequence in not a multiple of 3.
    if param.translate:
        import warnings
        from Bio import BiopythonWarning
        warnings.simplefilter('ignore', BiopythonWarning)

    # Interactive mode returns the origin.
    if param.genome or param.inter:
        recs = jsonrec.get_origin(item, param=param)
    elif param.protein:
        recs = jsonrec.get_translation_records(item, param=param)
    elif param.fasta:
        recs = jsonrec.get_feature_records(item, param=param)
    else:
        # Sanity check to hide tacit falltrough.
        raise Exception("this is not a valid option")

    return recs


def print_fasta(recs):
    """
    Prints fasta records
    """
    for rec in recs:
        print(rec.format("fasta"))


def fasta_view(params):
    """
    Converts data to fastya
    """

    for param in params:

        # Stop when data was not found.
        if not param.json:
            utils.error(f"data not found: {param.acc}")

        # Each data may have multiple entries.
        for item in param.json:

            # Get the fasta for each entry.
            recs = get_fasta(item, param=param)

            # Print the fasta records.
            print_fasta(recs)


@plac.pos("data", "data names")
@plac.flg('genome', "use the origin (genome) sequence", abbrev='g')
@plac.flg('fasta', "produce FASTA format", abbrev='F')
@plac.flg('protein', "operate on proteins", abbrev='p')
@plac.flg('translate', "translate DNA to protein", abbrev='T')
@plac.flg('transcribe', "transcribe DNA to RNA", abbrev='X')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type", abbrev="t")
@plac.opt('start', "start coordinate", abbrev="s")
@plac.opt('name', "select features by name", abbrev="n")
@plac.opt('id_', "select feature by id", abbrev="u")
@plac.opt('end', "end coordinate", abbrev="e")
@plac.opt('gene', "select features associated with a gene name", abbrev="G")
@plac.opt('match', "select features by rexep match")
@plac.flg('inter', "interactive (data from command line)", abbrev='i')
@plac.flg('reverse', "reverse sequence", abbrev='R')
@plac.flg('complement', "complement sequence", abbrev='C')
@plac.flg('revcomp', "reverse complement sequence", abbrev='r')
@plac.flg('verbose', "verbose mode")
def run(genome=False, fasta=False, protein=False, translate=False, transcribe=False, reverse=False,
        complement=False, revcomp=False, seqid='', start='', end='', type='', gene='', name='', match='', id_='',
        inter=False, verbose=False, *data):
    """
    Produces FASTA representations for data.
    """

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Reset counter (needed for consistency during testing).
    jsonrec.reset_counter()

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
                          acc=acc, translate=translate, reverse=reverse, uid=id_,
                          complement=complement, genome=genome, name=name, inter=inter,
                          fasta=fasta, type=type, gene=gene, regexp=match, transcribe=transcribe)

        # Fill the json data for the parameter if not an update
        p.json = storage.get_json(p.acc, seqid=seqid, inter=inter)
        return p

    # Each accession gets a parameter list.
    params = list(map(make_param, data))

    # Render the FASTA view.
    fasta_view(params)
