import sys
from itertools import repeat

from biorun import utils, models, parser
from biorun.libs import placlib as plac


@plac.pos("fnames")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.flg("diff", "output differences")
@plac.flg("table", "output in tabular format", abbrev="T")
@plac.flg("paired", "fasta input is pairwise", abbrev="p")
@plac.flg("vcf", "output vcf")
def run(start='', end='', diff=False, vcf=False, table=False, paired=False, *fnames):
    # fname = '../test/data/mafft.fa'

    # Parse the input
    recs = parser.get_records(fnames)

    recs = map(parser.make_seqrec, recs)

    alns = []
    try:

        if paired:
            # Sequences are target/query format
            stream = zip(recs, recs)
        else:
            # First sequence compared to all others
            targets = repeat(next(recs))
            stream = zip(targets, recs)

        for target, query in stream:

            # Sanity check.
            if len(target) != len(query):
                utils.error(
                    f"length of query and target do not match: {query.id}={len(query)}, {target.id}={len(target)}")

            # Parse start and end into user friendly numbers.
            if start or end:
                start = utils.parse_number(start)
                end = utils.parse_number(end)
                query.seq = query.seq[start:end]
                target.seq = target.seq[start:end]

            aln = models.Alignment(query=query, target=target,)

            alns.append(aln)


    except StopIteration as exc:
        utils.error(f'Input must have at least two FASTA sequences')
        sys.exit(1)

    if vcf:
        models.format_vcf(alns)
    elif diff:
        models.format_mutations(alns)
    elif table:
        models.format_table(alns)
    else:
        models.format_pairwise(alns)

if __name__ == '__main__':
    plac.call(run)
