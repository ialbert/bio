"""
Generates GFF outputs from a JSON record.
"""
from biorun import utils, const
from biorun.models import jsonrec
from itertools import count





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
            feats = jsonrec.filter_features(feats, start=param.start, end=param.end, gene=param.gene, ftype=param.type,
                                            regexp=param.regexp, name=param.name)

            # Generate the gff output
            for feat in feats:
                for values in feature2gff(feat, anchor=anchor, allow_parent=not(param.type)):
                    values = map(str, values)
                    print("\t".join(values))
