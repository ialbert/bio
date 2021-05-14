import gzip
import os
import sys
from collections import OrderedDict, defaultdict
from itertools import count
from uuid import uuid4

import biorun.libs.placlib as plac
from biorun import utils
from biorun.alias import ALIAS

logger = utils.logger

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError as exc:
    logger.error(f"{__name__}", stop=False)
    logger.error(f"{exc}", stop=False)
    logger.error(f"This software requires biopython.")
    logger.error(f"Try: conda install biopython>=1.78")
    sys.exit()

# GenBank terms to remap according to sequence ontology.
SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}

class Record:
    SOURCE = "source"
    def __init__(self, id, seqid=None, ftype=None, ann=None, meta=None, name=None, description=None):
        self.id = str(id)
        self.seqid = seqid or self.id
        self.name = name or ''
        self.description = description or ''
        self.type = ftype or self.SOURCE
        self.start = self.end = self.strand = 0
        self.seq_data = None
        self.seq_func = lambda: Seq('A')
        self.ann = ann or {}
        self.meta = meta or {}

    @property
    def seq(self):
        if self.seq_data is None:
            self.seq_data = self.seq_func()
        return self.seq_data

    def __str__(self):
        return f"Record: {self.id}, {self.type}"

def fasta(rec):
    seqrec = SeqRecord(seq=rec.seq, id=rec.id, description=rec.description)
    print(seqrec.format("fasta"), end='')

def first(data, key, default=""):
    # First element of a list value that is stored in a dictionary by a key.
    return data.get(key, [default])[0]

def make_uuid(ftype, size=6):
    """
    Unique id for a type.
    """
    uuid = str(uuid4())[:size]
    return f'{ftype}-{uuid}'

def find_name(ftype, ann):
    """
    Attempts to generate an unique id, a name and a description from annotations.
    """
    uid = desc = name = ''



    if ftype == 'gene':
        name = first(ann, "gene")
        desc = first(ann, "locus_tag")
    elif ftype == 'CDS':
        gene = first(ann, "gene")
        name = first(ann, "protein_id")
        desc = f"gene={gene} " + first(ann, "product")
    elif ftype == 'mRNA':
        name = first(ann, "transcript_id")
    elif ftype == "exon":
        name = first(ann, "gene")
        uid = make_uuid(ftype)
    else:
        desc = first(ann, "product")

    desc = f"{ftype} {desc}"
    name = name or make_uuid(ftype)
    uid = uid or name or make_uuid("uid")
    return uid, name, desc


def parse(fname):

    recs = SeqIO.parse(fname, format="genbank")

    for seqrec in recs:

        # Segrecord level annotations
        meta = dict(seqrec.annotations)

        seqid = ALIAS.get(seqrec.id, seqrec.id)

        for feat in seqrec.features:

            # Remap types to valid SO trams.
            ftype = SEQUENCE_ONTOLOGY.get(feat.type, feat.type)

            # Feature level annotations.
            ann = dict(feat.qualifiers)

            # Source features are special case
            if ftype == Record.SOURCE:
                uid, name, description = seqid, seqrec.name, seqrec.description
                func = lambda: seqrec.seq
            else:
                uid, name, description = find_name(ftype, ann=ann)
                func = lambda: feat.extract(seqrec.seq)

            # Remap the ids as well.
            uid = ALIAS.get(uid, uid)

            # Build the sequence record
            rec = Record(id=uid, seqid=seqid, ann=ann, meta=meta, name=name, description=description, ftype=ftype)
            rec.seq_func = func

            yield rec



def sequence_slicer(start=0, end=None):
    def func(rec):
        if start or end:
            seqlen = len(rec.seq)

            endx = seqlen if end is None else end
            endx = seqlen + endx if endx < 0 else endx

            # Zero based shift
            if start < 0:
                startx = seqlen + start + 1
            else:
                startx = start + 1

            rec.description = f"[{startx}:{endx}] {rec.description}"
            rec.seq_data = rec.seq[start:end]
        return rec

    return func


def type_selector(ftype):
    types = set(ftype.split(","))

    def func(rec):
        return rec.type in types if ftype else True

    return func


def gene_selector(name):
    genes = set(name.split(","))

    def func(rec):
        return first(rec.annot, "gene") in genes if name else True

    return func


def name_selector(name):
    names = set(name.split(","))

    def func(rec):
        return rec.id in names if name else True

    return func


def seqid_selector(seqid):
    targets = set(seqid.split(","))

    def func(rec):
        return rec.seqid in targets if seqid else True

    return func


def translate_recs(flag):
    def func(rec):
        if flag:
            endx = len(rec.obj.seq) // 3
            rec.obj.seq = rec.obj.seq[0:endx * 3].translate()
        return rec

    return func


def source_only(flag):

    def func(rec):
        return rec.type == Record.SOURCE if flag else True

    return func


def has_translation(rec):
    return "translation" in rec.annot


def get_translation(rec):
    rec.seq_data = Seq(rec.annot.get("translation")[0])
    return rec

def transform(recs, start=0, end=None, ftype=None, name=None, gene=None, protein=None, translate=None):

    recs = map(sequence_slicer(start=start, end=end), recs)

    return recs

if __name__ == '__main__':
    rec = Record(id=1)
    fasta(rec)

    fname = "../work/genomes.gb"
    recs = parse(fname)

    recs = transform(recs, end=10)

    for rec in recs:
        fasta(rec)

