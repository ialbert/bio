"""
Generates GFF outputs from a JSON record.
"""
from biorun import utils, const
from biorun.models import jsonrec
from itertools import count


def make_attr(feat, uid='', pid='', color=None):
    """
    Creates GFF style attributes from JSON fields.
    """

    # The feature name.
    name = feat['name']

    data = [ ]

    # Add ID and parent if these exists.
    if uid:
        data.append(f"ID={uid}")
    if pid:
        data.append(f"Parent={pid}")

    # Add the data name.
    data.append(f"Name={name}")

    # Other gff attributes.
    pairs = [(k, v) for k, v in feat.items() if k not in const.SKIP_GFF_ATTR]

    # Fill in all the fields.
    for key, value in pairs:
        data.append(f"{key}={value[0]}")

    # Attach a color to the feature.
    if color:
        data.append(f"color={color}")

    return ";".join(data)



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

    color = const.COLOR_FOR_TYPE.get(ftype)

    # Child type

    ctype = ftype

    pid = None

    # Handle hiearachical relationships.
    if ftype == 'mRNA' or ftype== 'CDS':
        pid = uid

        # Set the parent/child types.
        if ftype == 'mRNA':
            ptype, ctype = 'mRNA', 'exon'
        else:
            ptype, ctype = 'region', 'CDS'

        # Build the attributes for the parent.
        attr = make_attr(feat, color=color, uid=uid)

        # Generate the parent entry.
        data = [anchor, ".", ptype, feat['start'], feat['end'], ".", strand, phase, attr]

        # Separate mRNA track
        yield data

    # Process the individual locations.
    for start, end, strand in feat["location"]:
        strand = "+" if strand > 0 else "-"
        attr = make_attr(feat, color=color, pid=pid)
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
            collect = []
            for feat in feats:
                for values in feature2gff(feat, anchor=anchor):
                    values = map(str, values)
                    print("\t".join(values))
