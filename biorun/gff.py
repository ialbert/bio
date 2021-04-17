"""
Generates GFF outputs from a JSON record.
"""
import sys, json
import biorun.libs.placlib as plac
from biorun import utils, jsonx
from urllib.parse import  quote

# Module level logger.
logger = utils.logger

# Associates a color to a feature type.
COLOR_FOR_TYPE = {
    "five_prime_UTR": "#cc0e74",
    "three_prime_UTR": "#cc0e74",
    "stem_loop": "#fa7f72",
    "mature_protein_region": "#CBAEBB",
    "region": "#CECECE",
    "mRNA": "#799351",
    "gene": "#cb7a77",
    "transcript": "#79a3b1",
    "tRNA": "#a685e2",
    "ncRNA": "#fca3cc",
    "mobile_element": "#efd9d1",
    "mRNA_region":"#7a77cb",
}
#
# The GFF attributes generated for a source type.
#
SOURCE_ATTRIBUTES = [
    "mol_type", "isolate", "db_xref", "organism", "country", "collection_date"
]

# GFF attributes filled for each feature other than "source"
GFF_ATTRIBUTES = [
    "gene", "protein_id", "product", "db_xref", "function",
]

# Attributes that should not be added
SKIP_GFF_ATTR = {"id", "parent_id", "name", "type", "start", "end", "location", "translation", "strand", "operator"}

def make_attr(feat, color=None, pid=None):
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

    # Add all attributes to tGFF.
    pairs = [(k, v) for k, v in feat.items() if k not in SKIP_GFF_ATTR]

    # Fill in all the fields.
    for key, value in pairs:
        val = quote(str(value[0]))
        data.append(f"{key}={val}")

    # Attach a color to the feature.
    if color:
        data.append(f"color={color}")

    text = ";".join(data)

    return  text



def feature2gff(feat, seqid):
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
    color = COLOR_FOR_TYPE.get(ftype)

    # Make the attributes
    attr = make_attr(feat, color=color)

    # Create the GFF record.
    data = [ seqid, ".", ftype, start, end, ".", gffstrand, phase, attr]

    yield data

@plac.flg("features", "convert the features")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "filter for a sequence id")
@plac.opt("name", "filter for a sequence name")
@plac.opt("gene", "filter for a gene name")
@plac.flg("proteins", "extract embedded protein sequences")
@plac.flg("translate", "translate DNA sequences", abbrev='R')
def run(features=False, proteins=False, translate=False,
        start='0', end='0', type_='', id_='', name='', gene=''):
    """
    Converts data to fastyaq
    """

    # Load the data
    data = json.loads(sys.stdin.read())

    # Parse start and end into user friendly numbers.
    start = utils.parse_number(start)
    end = utils.parse_number(end)
    ftype = type_
    seqid = id_

    feats = jsonx.select_features(data)

    print("##gff-version 3")
    for seqid, feat in feats:
        for values in feature2gff(feat, seqid=seqid):
            values = map(str, values)
            print("\t".join(values))

