import sys

from biorun import utils, models, parser
from biorun.libs import placlib as plac


@plac.pos("fnames")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.flg("diff", "output differences")
@plac.flg("vcf", "output vcf")
def run(start='', end='', diff=False, vcf=False, *fnames):
    # fname = '../test/data/mafft.fa'

    # Parse the input
    recs = parser.get_records(fnames)

    recs = iter(recs)

    try:
        query = next(recs)

        for target in recs:

            # Sanity check.
            if len(target) != len(query):
                utils.error(f"# length of query and target do not match: {len(query)}, {len(target)}")

            # Parse start and end into user friendly numbers.
            if start or end:
                start = utils.parse_number(start)
                end = utils.parse_number(end)
                query.seq = query.seq[start:end]
                target.seq = target.seq[start:end]

            aln = models.Alignment(query=query, target=target)

            alns = [aln]

            if vcf:
                models.format_vcf(alns)
            elif diff:
                models.format_diffs(alns)
            else:
                models.format_pairwise(alns)

    except StopIteration as exc:
        utils.error(f'Input must have at least two FASTA sequences')
        sys.exit(1)





if __name__ == '__main__':
    plac.call(run)
