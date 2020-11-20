"""
Generates GFF outputs from a JSON record.
"""
from biorun import utils, const
from biorun.models import jsonrec


def feature2gff(feat, anchor):
    """
    Returns a SeqRecord as an 11 element  GFF3 list .
    """

    ftype = feat['type']


    # TODO: is this the phase?
    #phase = feat.get("codon_start", [1])[0] - 1
    phase = "."

    for start, end, strand in feat["location"]:
        strand = feat['strand']
        attr = jsonrec.make_attr(feat)

        # Generate a color for the feature type.
        color = const.COLOR_FOR_TYPE.get(ftype)
        if color:
            attr += f";color={color}"

        strand = "+" if strand else "-"

        data = [anchor, ".", ftype, start, end, ".", strand, phase, attr]

        yield data


def gff_view(params):
    """
    Converts data to fastya
    """

    print("##gff-version 3")

    for param in params:

        # Stop when data was not found.
        if not param.json:
            utils.error(f"data not found: {param.name}")

        # Each data may have multiple entries.
        for item in param.json:

            # Pull out the features.
            feats = item[const.FEATURES]

            # The name of the GFF anchor.
            anchor = param.seqid or item['id']

            # Subselect by coordinates.
            feats = jsonrec.filter_features(feats, start=param.start, end=param.end, gene=param.gene, ftype=param.type,
                                            regexp=param.regexp)

            # Generate the gff output
            for feat in feats:
                for values in feature2gff(feat, anchor=anchor):
                    values = map(str, values)
                    print("\t".join(values))
