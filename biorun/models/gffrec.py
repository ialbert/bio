"""
Generates GFF outputs from a JSON record.
"""
import biorun.libs.placlib as plac
from biorun.models import jsonrec
from biorun import utils, const, fetch, objects

# Module level logger.
logger = utils.logger


def make_attr(feat, color=None):
    """
    Creates GFF style attributes from JSON fields.
    """

    # The feature name.
    name = feat['name']
    uid = feat['id']
    pid = feat.get("parent_id")

    data = [f"ID={uid}"]

    # Add parent id
    if pid:
        data.append(f"Parent={pid}")

    # Add the data name.
    data.append(f"Name={name}")

    # Fill in all the fields.
    for key, value in pairs:
        data.append(f"{key}={value[0]}")

    # Attach a color to the feature.
    if color:
        data.append(f"color={color}")

    return ";".join(data)



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

        # GFF is a interval (record mode).
        param.record = True

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

