"""
Generates GFF outputs from a JSON record.
"""
from biorun import utils, const
from biorun.models import jsonrec

def get_color(ftype):
    """
    Generates a color for a type.
    """
    color = const.COLOR_FOR_TYPE.get(ftype)
    return f";color={color}" if color else ''

def feature2gff(feat, anchor):
    """
    Returns a SeqRecord as an 11 element  GFF3 list .
    """

    uid = feat['id']
    ftype = feat['type']
    strand = feat['strand']

    # Reformat the strand
    strand = "+" if strand > 0 else "-"

    # TODO: is this the phase?
    #phase = feat.get("codon_start", [1])[0] - 1
    phase = "."

    # Produce a GFF representation for the original feature.
    attr1 = jsonrec.make_attr(feat)

    # Color the attribute.
    attr1 = attr1 + get_color(ftype)

    # Create the parent attribute.
    data = [anchor, ".", ftype, feat['start'], feat['end'], ".", strand, phase, attr1]

    yield data

    # Iterate over locations for mRNA to generate exons
    if ftype == 'mRNA':

        for start, end, strand in feat["location"]:

            strand = "+" if strand > 0 else "-"

            ctype = 'exon'

            attr2 = jsonrec.make_attr(feat)

            attr2 = attr2 + get_color(ctype)

            data = [anchor, ".", ctype, start, end, ".", strand, phase, attr2]

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
