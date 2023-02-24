from biorun.libs import placlib as plac
import pandas as pd
from gprofiler import GProfiler
import pickle

@plac.opt('counts', "input counts", abbrev='c')
@plac.opt('organism', "input counts", abbrev='d')
@plac.opt('colname', "gene id column name", abbrev='n')
@plac.opt('pval_cutoff', "pvalue cutoff", abbrev='t', type=float)
@plac.opt('pval_column', "pvalue column name", abbrev='p')
@plac.opt('output', "pvalue column name", abbrev='o')
def run(counts="edger.csv", organism='mmusculus', colname='gene', pval_cutoff=0.05, pval_column='FDR', output='gprofiler.csv'):
    """
    Runs the g:Profiler tool on a csv file where one column contains gene names and some column contains pvalues.

    Filters the p values by a threshold, then submits the gene names to g:GOSt.

    Valid organisms: https://biit.cs.ut.ee/gprofiler/page/organism-list
    """
    ct = pd.read_csv(counts)

    ct = ct[ct[pval_column] < pval_cutoff]

    query = ct[colname].tolist()

    def strip_dot(x):
        return x.split('.')[0]

    # Get rid of version numbers if these exists
    query = list(map(strip_dot, query))

    print(f"# Running g:Profiler")
    #print(f"# https://biit.cs.ut.ee/gprofiler/gost")
    print(f"# Counts: {counts}")
    print(f"# Organism: {organism}")
    print(f"# Name column: {colname}")
    print(f"# Pval column: {pval_column} < {pval_cutoff}")
    print(f"# Gene count: {len(query)}")
    print(f"# Genes: {','.join(query[:5])},[...]")
    print(f"# Submitting to gProfiler")

    # Submit the query to g:Profiler.
    gp = GProfiler(return_dataframe=True)
    res = gp.profile(organism=organism, query= query)

    res = res.sort_values('source')
    res.to_csv(output, index=False)
    #print(res)
    print(f"# Found {len(res)} functions")
    print(f"# Output: {output}")

if __name__ == '__main__':
    plac.call(run)
