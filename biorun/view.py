"""
Handles functionality related to data storege.
"""
from biorun import const, utils, objects, fetch
from biorun.models import jsonrec, fastarec, gffrec
import biorun.libs.placlib as plac

# Module level logger.
logger = utils.logger


@plac.pos("data", "data names")
@plac.flg('protein', "operate on proteins", abbrev='p')
@plac.flg('translate', "translate DNA to protein", abbrev='T')
@plac.flg('transcribe', "transcribe DNA to RNA", abbrev='X')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type", abbrev="t")
@plac.opt('start', "start coordinate", abbrev="s")
@plac.opt('name', "select features by name", abbrev="n")
@plac.opt('id_', "select feature by id", abbrev="u")
@plac.opt('end', "end coordinate", abbrev="e")
@plac.opt('gene', "select features associated with a gene name", abbrev="g")
@plac.opt('match', "select features by rexep match")
@plac.flg('inter', "interactive (data from command line)", abbrev='i')
@plac.flg('reverse', "reverse sequence", abbrev='R')
@plac.flg('complement', "complement sequence", abbrev='C')
@plac.flg('revcomp', "reverse complement sequence", abbrev='r')
@plac.flg('features', "operate on features", abbrev='f')
@plac.flg('fasta', "produce FASTA format", abbrev='F')
@plac.flg('gff', "produce GFF format", abbrev='G')
@plac.flg('json', "produce JSON format", abbrev='J')
@plac.flg('genbank', "produce GenBank format", abbrev='B')
@plac.flg('verbose', "verbose mode")
def run(protein=False, translate=False, transcribe=False, reverse=False,
        complement=False, revcomp=False, seqid='', start='', end='', type='', gene='', name='', match='', id_='',
        inter=False, features=False, fasta=False,  gff=False, json=False, genbank=False, verbose=False, *data):
    """
    Produces FASTA representations for data.
    """

    # Turn on features if some parameters are present.
    features = features or (type or name or match or id_ or protein)

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Reset counter (needed for consistency during testing).
    jsonrec.reset_counter()

    # Check that data have no dashes.
    utils.no_dash(data)

    def make_param(acc):
        """
        Creates a parameter for each accession.

        """
        # Set the verbosity
        utils.set_verbosity(logger, level=int(verbose))

        # A simple wrapper class to carry all parameters around.
        p = objects.Param(start=start, end=end, seqid=seqid, protein=protein, revcomp=revcomp,
                          acc=acc, translate=translate, reverse=reverse, uid=id_, gff=gff,
                          complement=complement, name=name, inter=inter, features=features,
                          fasta=fasta, type=type, gene=gene, regexp=match, transcribe=transcribe)

        # Fill the json data for the parameter if not an update
        p.json = fetch.get_json(p.acc, seqid=seqid, inter=inter)

        return p

    params = list(map(make_param, data))

    if fasta:
        fastarec.fasta_view(params)
    elif gff:
        gffrec.gff_view(params)
    elif json:
        jsonrec.json_view(params)
    elif genbank:
        fetch.genbank_view(params)
    else:
        fastarec.fasta_view(params)