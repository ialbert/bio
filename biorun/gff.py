import plac

from biorun import convert


@plac.pos("data", "input data")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "filter for a sequence id")
@plac.opt("name", "filter for a sequence name")
@plac.opt("gene", "filter for a gene name", abbrev='g')
@plac.opt("alias", "remap sequence ids")
@plac.pos("fnames", "input files")
def run(start='1', end='', type_='', id_='', name='', gene='', alias=None, *fnames):
    convert.run(start=start, end=end, type_=type_, id_=id_, name=name, gene=gene, alias=alias, fasta=False, fnames=fnames)


