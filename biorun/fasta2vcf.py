import plac
from itertools import *
from biorun import utils
from biorun import align


@plac.pos("fnames")
def run(*fnames):

    #fname = '../test/data/mafft.fa'
    streams = utils.open_streams(fnames=fnames)

    for stream in streams:

        try:
            recs = utils.fasta_parser(stream)
            ref = next(recs)
            query = next(recs)
        except StopIteration as exc:
            continue

        vcfdict = align.find_variants(ref, query)

        print('##fileformat=VCFv4.2')
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        print('##FILTER=<ID=PASS,Description="All filters passed">')
        print('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">')
        print(f'##contig=<ID={ref.name},length={len(ref.seq.strip("-"))},assembly={ref.name}>')
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{query.name}")

        for value in vcfdict.values():
            print("\t".join(value))

if __name__ == '__main__':
    plac.call(run)
