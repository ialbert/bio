import sys
from biorun.libs import placlib as plac
from biorun import convert


@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "exact match on a sequence id")
@plac.opt("match", "regexp match on a sequence id")
@plac.opt("gene", "filter for a gene name", abbrev='g')
@plac.opt("fields", "table fields (default: id,size)", abbrev='f')
@plac.opt("rename", "rename sequence ids")
@plac.opt("olap", "overlap with coordinate")
@plac.pos("fnames", "input files")
def run(start='1', end='', type_='', id_='', match='', gene='',
        rename='', olap='', fields="id,type,size",  *fnames):
    """
    Generates tabular output from data.
    """


    convert.run(protein=False, translate=False, start=start, table=True, fields=fields,
                end=end, type_=type_, id_=id_, revcomp=False, features=None, olap=olap,
                match=match, gene=gene, rename=rename, fasta=True, fnames=fnames)

