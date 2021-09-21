from biorun.libs import placlib as plac
import sys
from itertools import *
from biorun import utils, models, parser


@plac.pos("fnames")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
def run(start='', end='', *fnames):


    #fname = '../test/data/mafft.fa'

    # Parse the input
    recs = parser.get_records(fnames)

    recs = iter(recs)

    try:
        query = next(recs)
        target = next(recs)

        # Parse start and end into user friendly numbers.
        if start or end:
            start = utils.parse_number(start)
            end = utils.parse_number(end)
            query.seq = query.seq[start:end]
            target.seq = target.seq[start:end]

    except StopIteration as exc:
        utils.error(f'Input must have at least two FASTA sequences')
        sys.exit(1)

    vcfdict = models.find_variants(query, target)

    print('##fileformat=VCFv4.2')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FILTER=<ID=PASS,Description="All filters passed">')
    print('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">')
    print(f'##contig=<ID={target.name},length={len(target.seq.strip("-"))},assembly={target.name}>')
    print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{query.name}")

    for value in vcfdict.values():
        print("\t".join(value))

    if len(target)!=len(query):
        utils.error(f"# length of query and target do not match: {len(query)}, {len(target)}")

if __name__ == '__main__':
    plac.call(run)
