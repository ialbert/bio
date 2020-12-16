
from goatools.cli.find_enrichment import GoeaCliArgs
from goatools.cli.find_enrichment import GoeaCliFnc

from biorun import utils
from biorun.libs import placlib as plac

logger = utils.logger


def build_population():
    return


@plac.pos('filenames', 'data/study data/population data/association')
@plac.opt('--study', type=str, 'Gene names in a study')
@plac.opt('--population', type=str, 'Gene names in population (or other study if --compare is specified)')
@plac.opt('--association', type=str, 'An association file that maps a gene name to a GO category.')
@plac.opt('--taxid', type=int, "When using NCBI's gene2go annotation file, specify desired taxid")
@plac.opt('--alpha', type=float, help='Test-wise alpha for multiple testing')
@plac.opt('--pval', type=float, help="""Only print results with uncorrected p-value < PVAL.  
                                        Print all results, significant and otherwise, by setting --pval=1.0""")
@plac.opt('--pval_field', type=str,
               help='Only print results when PVAL_FIELD < PVAL.')
@plac.opt('--outfile', type=str,
               help='Write enrichment results into xlsx or tsv file')
@plac.opt('--ns', type=str,
               help='Limit GOEA to specified branch categories. '
                    'BP=Biological Process; '
                    'MF=Molecular Function; '
                    'CC=Cellular Component')
@plac.opt('--id2sym', type=str,
               help='ASCII file containing one geneid and its symbol per line')
@plac.opt('--outfile_detail', type=str,
               help=('Write enrichment results into a text file \n'
                     'containing the following information: \n'
                     '1) GOEA GO terms, grouped into sections \n\n'
                     '2) List of genes and ASCII art showing section membership \n'
                     '3) Detailed list of each gene and GO terms w/their P-values \n'))
@plac.flg('--compare',
               help="the population file as a comparison group. if this "
                    "flag is specified, the population is used as the study "
                    "plus the `population/comparison`")
@plac.opt('--ratio', type=float,
               help="only show values where the difference between study "
                    "and population ratios is greater than this. useful for "
                    "excluding GO categories with small differences, but "
                    "containing large numbers of genes. should be a value "
                    "between 1 and 2. ")
@plac.flg('--indent', help="indent GO terms")
@plac.opt('--obo', type=str,
               help="Specifies location and name of the obo file")
@plac.flg('--no_propagate_counts',
               help="Do not propagate counts to parent terms")
# no -r:   args.relationship == False
# -r seen: args.relationship == True
@plac.flg('--relationship',
               help="Propagate counts up all relationships,")
def run(*taxids):

    # Set the verbosity
    utils.set_verbosity(logger, level=int(0))

    # Load study, population, associations, and GoDag. Run GOEA.
    obj = GoeaCliFnc()

    # Reduce results to significant results (pval<value)
    results_specified = obj.get_results()

    # Print results in a flat list
    obj.prt_results(results_specified)

    pass

