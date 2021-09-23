import sys
from biorun.libs import placlib as plac
from biorun import convert


@plac.pos("data", "input data")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "filter for a sequence id")
@plac.opt("name", "filter for a sequence name")
@plac.opt("gene", "filter for a gene name", abbrev='g')
@plac.flg("protein", "operate on the protein sequences", abbrev='p')
@plac.flg("translate", "translate DNA ", abbrev='T')
@plac.flg("revcomp", "reverse complement DNA", abbrev='R')
@plac.opt("alias", "remap sequence ids")
@plac.opt("frame", "reading frame", type=int, choices=[1,2,3,-1,-2,-3], abbrev='F')
@plac.pos("fnames", "input files")
def run(start='1', end='', type_='', id_='', name='', gene='',
        alias='', protein=False, translate=False, revcomp=False, frame=1, *fnames):

    if frame > 1:
        start = frame

    if frame < 1:
        start = -frame
        revcomp = True

    convert.run(protein=protein, translate=translate, start=start,
                end=end, type_=type_, id_=id_, revcomp=revcomp,
                name=name, gene=gene, alias=alias, fasta=True, fnames=fnames)

