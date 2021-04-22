"""
Converts across formats
"""
import sys, json, os, pathlib, gzip
import biorun.libs.placlib as plac
from biorun import utils
from biorun.alias import ALIAS
from functools import partial

# Module level logger.
logger = utils.logger

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import Reference, CompoundLocation, FeatureLocation
except ImportError as exc:
    logger.error(f"{__name__}", stop=False)
    logger.error(f"{exc}", stop=False)
    logger.error(f"This software requires biopython.")
    logger.error(f"Try: conda install biopython>=1.78")
    sys.exit()

def is_fasta(fname):
    exts = "fa fasta fa.gz fasta.gz".split()
    found = filter(lambda x:fname.endswith(x), exts)
    return any(found)

def is_genbank(fname):
    exts = "gb gb.gz genbank genbank.gz gpff gpff.gz".split()
    found = filter(lambda x:fname.endswith(x), exts)
    return any(found)

def parse(fname):
    """
    Parses a filename with the appropriate readers.
    """
    stream = gzip.open(fname) if fname.endswith("gz") else open(fname)
    if is_fasta(fname):
        recs = SeqIO.parse(stream, format="fasta")
    elif is_genbank(fname):
        recs = SeqIO.parse(stream, format="genbank")
    else:
        logger.error(f"file extension not recognized: {fname}")
        sys.exit()

    return recs

def make_record(text, seqid=1, locus="", desc=""):
    seq = Seq(text)
    rec = SeqRecord(seq=seq, id=seqid, name=locus, description=desc)
    return rec

def read_input(fname, store=None, interactive=False):
    """
    Attempts load the correct input.
    """

    # Item is a valid file.
    if os.path.isfile(fname):
        data = parse(fname)
        return data

    # Generate an interactive data.
    if interactive:
        data = make_record(text=fname)
        return data

    # Invalid data.
    logger.error(f"file not found: {fname}")
    sys.exit()


@plac.pos("data", "input data")
@plac.flg("features", "convert the features", abbrev='F')
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
def run(features=False, proteins=False, translate=False, gff_=False, fasta_=False,
        start='1', end=None, type_='', id_='', name='', gene='', *fnames):
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
    fasta_ = False if (gff_ and not fasta_) else True

    def fasta_formatter(rec, start=None, end=None):
        print(rec.format("fasta"), end='')

    def remapper(rec):
        rec.id = ALIAS.get(rec.id, rec.id)
        return rec

    def slicer(rec):
        if start or end:
            startx = start - 1
            endx = end if end else len(rec.seq)
            rec.description = f"{rec.description} [{start+1}:{endx}]"
            rec.seq = rec.seq[startx:end]
        return rec

    formatter = fasta_formatter

    for fname in fnames:
        recs = read_input(fname, interactive=False)
        recs = map(remapper, recs)
        recs = map(slicer, recs)
        for rec in recs:
            formatter(rec)

    sys.exit()

    '''
    for datreca in inputs:
        if fasta_:
            recs = jsonx.select_records(data, features=features, proteins=proteins, translate=translate,
                                        start=start, end=end, ftype=ftype, seqid=seqid, name=name, gene=gene)
            for rec in recs:


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
    '''