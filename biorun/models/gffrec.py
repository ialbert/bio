"""
Generates GFF outputs from a JSON record.
"""
import biorun.libs.placlib as plac
from biorun.models import jsonrec
from biorun import utils, const, storage, objects

# Module level logger.
logger = utils.logger

def feature2gff(feat, anchor, allow_parent=True):
    """
    Returns a SeqRecord as an 11 element  GFF3 list .
    """

    ftype = feat['type']
    strand = feat['strand']
    start = feat['start']
    end = feat['end']

    # Reformat the strand
    gffstrand = "+" if strand > 0 else "-"

    # TODO: is this the phase?
    #phase = feat.get("codon_start", [1])[0] - 1
    phase = "."

    # The color for the feature.
    color = const.COLOR_FOR_TYPE.get(ftype)

    # Make the attributes
    attr = jsonrec.make_attr(feat, color=color)

    # Create the GFF record.
    data = [anchor, ".", ftype, start, end, ".", gffstrand, phase, attr]

    yield data


def gff_view(params):
    """
    Converts data to fastya
    """

    print("##gff-version 3")

    for param in params:

        # Stop when data was not found.
        if not param.json:
            utils.error(f"data not found: {param.acc}")

        # Each data may have multiple entries.
        for item in param.json:

            # Pull out the features.
            feats = jsonrec.get_json_features(item)

            # The name of the GFF anchor.
            anchor = param.seqid or item['id']

            # Subselect by coordinates.
            feats = jsonrec.filter_features(feats, param=param)

            # Generate the gff output
            for feat in feats:
                for values in feature2gff(feat, anchor=anchor, allow_parent=not(param.type)):
                    values = map(str, values)
                    print("\t".join(values))


@plac.pos("data", "data names")
@plac.flg('protein', "operate on proteins", abbrev='p')
@plac.opt('seqid', "set the sequence id", abbrev='S')
@plac.opt('type', "select feature by type", abbrev="t")
@plac.opt('start', "start coordinate", abbrev="s")
@plac.opt('name', "select feature by name", abbrev="n")
@plac.opt('id_', "select feature by id", abbrev="u")
@plac.opt('end', "end coordinate", abbrev="e")
@plac.opt('gene', "select features associated with a gene name", abbrev="G")
@plac.opt('match', "select features by rexep match")
@plac.flg('inter', "interactive (data from command line)", abbrev='i')
@plac.flg('gff', "activate gff subcommand", abbrev='g')
@plac.flg('verbose', "verbose mode")
def run(protein=False, seqid='', start='', end='', type='', gene='', name='', match='', id_='',
        inter=False, gff=False, verbose=False,  *data):
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

        # A simple wrapper class to carry all parameters around.
        p = objects.Param(acc=acc, start=start, end=end, seqid=seqid, protein=protein, name=name, inter=inter,
                          uid=id_, type=type, gene=gene, regexp=match)

        # Turn off sequence mode.
        p.genome = p.fasta = False

        # Fill the json data for the parameter if not an update
        p.json = storage.get_json(p.acc, seqid=seqid, inter=inter)
        return p

    # Each accession gets a parameter list.
    params = list(map(make_param, data))

    # Render the FASTA view.
    gff_view(params)