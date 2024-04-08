from biorun.libs import placlib as plac

from biorun import convert


@plac.pos("data", "input data")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "filter for a sequence id")
@plac.opt("match", "regexp match on a name")
@plac.opt("olap", "overlap with coordinate")
@plac.opt("gene", "filter for a gene name", abbrev='g')
@plac.opt("rename", "rename sequence ids")
@plac.pos("fnames", "input files")
def run(start='1', end='', type_='', id_='', match='', gene='', olap='', rename='', *fnames):
    convert.run(start=start, end=end, type_=type_, match=match, id_=id_, gene=gene, rename=rename, olap=olap, fasta=False, fnames=fnames)


