"""
Generates GFF outputs from a JSON record.
"""
from biorun import utils
from urllib.parse import  quote
from pprint import pprint
from itertools import count


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
SKIP_GFF_ATTR = {"id", "parent_id", "name", "type", "start", "comment", "references", "structured_comment",
                 "end", "location", "translation", "strand", "operator"}


def feature2gff(seqid, ftype, start, end, strand, uid, name, pid=None):
    """
    Returns a Record as an 11 element  GFF3 list .
    """
    # Reformat the strand
    strand = "+" if strand > 0 else "-"

    # TODO: is this the phase?
    #phase = feat.get("codon_start", [1])[0] - 1
    phase = "."

    # The color for the feature.
    color = COLOR_FOR_TYPE.get(ftype)

    # Attribute data
    attr = [ f"ID={uid}", f"Name={name}" ]
    if pid:
        attr.append( f"Parent={pid}")
    if color:
        attr.append(f"color={color}")

    # Build the attribute string
    attr = ";".join(attr)

    # Create the GFF record.
    data = [ seqid, ".", ftype, start, end, ".", strand, phase, attr]

    return data



counter = count(1)

def gff_formatter(rec):

    # Parent feature
    data = feature2gff(start=rec.start, end=rec.end, ftype=rec.type, uid=rec.id, name=rec.name, strand=rec.strand, seqid=rec.seqid, pid=None)
    line = "\t".join(map(str, data))

    # Parent id.
    pid = rec.id

    if rec.type == "mRNA":
        print (line)
        ftype = "exon"
    else:
        ftype = rec.type

    # Generate the locations
    for start, end, strand in rec.locations:
        name = rec.name
        uid = next(counter)

        data = feature2gff(start=start, end=end, ftype=ftype, uid=uid, name=name,
                           strand=strand, seqid=rec.seqid, pid=rec.id)
        line = "\t".join(map(str, data))
        print(line)

