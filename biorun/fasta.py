import sys
from biorun.libs import placlib as plac
from biorun import convert


@plac.pos("data", "input data")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "exact match on a sequence id")
@plac.opt("match", "regexp match on a sequence id")
@plac.opt("gene", "filter for a gene name", abbrev='g')
@plac.flg("protein", "operate on the protein sequences", abbrev='p')
@plac.flg("translate", "translate DNA ", abbrev='T')
@plac.flg("features", "extract the fasta for the features", abbrev='G')
@plac.flg("revcomp", "reverse complement DNA", abbrev='R')
@plac.opt("trim", "trim polyA tails (and leading/trailing Ns)", abbrev='A')
@plac.opt("rename", "rename sequence ids")
@plac.opt("olap", "overlap with coordinate")
@plac.opt("frame", "reading frame", type=int, choices=[1,2,3,-1,-2,-3], abbrev='F')
@plac.pos("fnames", "input files")
def run(start='1', end='', type_='', id_='', match='', gene='',
        rename='', protein=False, translate=False, revcomp=False, features=False, trim='', olap='', frame=1, *fnames):

    if frame > 1:
        start = frame

    if frame < 1:
        start = -frame
        revcomp = True

    convert.run(protein=protein, translate=translate, start=start,
                end=end, type_=type_, id_=id_, revcomp=revcomp, features=features, olap=olap, trim=trim,
                match=match, gene=gene, rename=rename, fasta=True, fnames=fnames)

