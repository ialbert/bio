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
@plac.flg("genome", "extract genome only ", abbrev='G')
@plac.flg("revcomp", "reverse complement DNA", abbrev='R')
@plac.opt("rename", "rename sequence ids")
@plac.opt("olap", "overlap with coordinate")
@plac.flg("size", "print accession and size only", abbrev='S')
@plac.opt("frame", "reading frame", type=int, choices=[1,2,3,-1,-2,-3], abbrev='F')
@plac.pos("fnames", "input files")
def run(start='1', end='', type_='', id_='', match='', gene='',
        rename='', protein=False, translate=False, revcomp=False, genome=False, olap='', size=False, frame=1, *fnames):

    if frame > 1:
        start = frame

    if frame < 1:
        start = -frame
        revcomp = True

    convert.run(protein=protein, translate=translate, start=start, size=size,
                end=end, type_=type_, id_=id_, revcomp=revcomp, genome=genome, olap=olap,
                match=match, gene=gene, rename=rename, fasta=True, fnames=fnames)

