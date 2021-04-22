"""
Converts across formats
"""
import sys, json, os, pathlib, gzip
import biorun.libs.placlib as plac
from biorun import utils, jsonx, reclib, gff
from functools import partial

# Module level logger.
logger = utils.logger

def is_type(fname, ext=[]):
    path = pathlib.Path(fname)
    expect = set(x.lower() for x in ext)
    found = set(x.lower().strip(".") for x in path.suffixes)
    if (found & expect):
        return True

    return False


def read_file(fname):
    """
    Returns a JSON representation from a file.
    """
    stream = gzip.open(fname) if fname.endswith("gz") else open(fname)
    if is_type(fname, ext=["json"]):
        data = json.load(stream)
    elif is_type(fname, ext=["fa", "fasta"]):
        data = jsonx.parse_stream(stream, type="fasta")
    elif is_type(fname, ext=["gb", "genbank", "gpff"]):
        data = jsonx.parse_stream(stream, type="genbank")
    else:
        logger.error(f"file extension not recognized: {fname}")
        sys.exit()
    return data


def read_inputs(fname, store=None, interactive=False):
    """
    Attempts load the correct input.

    Priorities:

    1. file loaded first
    2. then database


    """

    # Item is a valid file.
    if os.path.isfile(fname):
        data = read_file(fname)
        return data

    # Try to load a database
    if store:
        items = [item for item in store if item['id'] == fname]
        if items:
            return items

    # Generate an interactive data.
    if interactive:
        data = jsonx.make_jsonrec(fname)
        return data

    # Invalid data.
    logger.error(f"file not found: {fname}")
    sys.exit()


@plac.pos("data", "input data")
@plac.flg("features", "convert the features", abbrev='F')
@plac.flg("json_", "convert to json")
@plac.flg("fasta_", "convert to fasta")
@plac.flg("gff_", "convert to gff")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "filter for a sequence id")
@plac.opt("name", "filter for a sequence name")
@plac.opt("gene", "filter for a gene name", abbrev='G')
@plac.flg("proteins", "operate on the protein sequences", abbrev='P')
@plac.flg("translate", "translate DNA sequences", abbrev='R')
def run(features=False, proteins=False, translate=False, json_=False, gff_=False, fasta_=False,
        start='0', end='0', type_='', id_='', name='', gene='', *fnames):
    """
    Convert data to various formats
    """

    # Parse start and end into user friendly numbers.
    start = utils.parse_number(start)
    end = utils.parse_number(end)
    ftype = type_
    seqid = id_

    elems = seqid.split(":")
    if len(elems) == 2:
        seqid, gene = elems
        ftype = "CDS"


    # Default format is fasta if nothing is specified.
    fasta_ = False if ((json_ or gff_) and not fasta_) else True

    # File remapper.
    reader = partial(read_inputs, interactive=False)

    # Generate the input data.
    inputs = map(reader, fnames)

    for data in inputs:
        if fasta_:
            recs = jsonx.select_records(data, features=features, proteins=proteins, translate=translate,
                                        start=start, end=end, ftype=ftype, seqid=seqid, name=name, gene=gene)
            for rec in recs:
                print(rec.format("fasta"), end='')

        elif gff_:
            feats = jsonx.select_features(data, start=start, end=end, ftype=ftype, seqid=seqid, name=name, gene=gene)

            print("##gff-version 3")
            for seqid, feat in feats:
                for values in gff.feature2gff(feat, seqid=seqid):
                    values = map(str, values)
                    print("\t".join(values))

        elif json_:
            text = json.dumps(list(data), indent=4)
            print(text)
