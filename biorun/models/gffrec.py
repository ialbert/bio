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

    uid = feat['id']
    ftype = feat['type']
    strand = feat['strand']

    # Reformat the strand
    strand = "+" if strand > 0 else "-"

    # TODO: is this the phase?
    #phase = feat.get("codon_start", [1])[0] - 1
    phase = "."

    # The color for the feature.
    color = const.COLOR_FOR_TYPE.get(ftype)

    # Used only during hierachical data.
    ctype, pid = ftype, None

    # Handle hiearachical relationships.
    if ftype == 'mRNA' or ftype== 'CDS':

        # Child nodes track parent.
        pid = uid

        # Set the parent/child types.
        if ftype == 'mRNA':
            ptype, ctype = 'mRNA', 'exon'
        else:
            ptype, ctype = 'region', 'CDS'

        # Build the attributes for the parent.
        attr = jsonrec.make_attr(feat, color=color, uid=uid)

        # Generate the parent entry.
        data = [anchor, ".", ptype, feat['start'], feat['end'], ".", strand, phase, attr]

        # Output the parent track if allowed.
        if allow_parent:
            yield data

    # Generate an interval for each location.
    for start, end, strand in feat["location"]:
        strand = "+" if strand > 0 else "-"
        attr = jsonrec.make_attr(feat, color=color, pid=pid)
        data = [anchor, ".", ctype, start, end, ".", strand, phase, attr]

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
                for values in feature2gff(feat, anchor=anchor, allow_parent=not(param.type)):
                    values = map(str, values)
                    print("\t".join(values))
