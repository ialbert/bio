
from goatools.cli.find_enrichment import GoeaCliArgs
from goatools.cli.find_enrichment import GoeaCliFnc

from biorun import utils
from biorun.libs import placlib as plac

logger = utils.logger


@plac.pos('taxids', "verbose mode, prints more messages")
@plac.pos('filenames','data/study data/population data/association')
@plac.opt('--annofmt',type=str, choices=['gene2go', 'gaf', 'gpad', 'id2gos'],
          'Annotation file format. Not needed if type can be determined using filename')
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
@plac.opt('--sections',  type=str,
               help=('Use sections file for printing grouped GOEA results. '
                     'Example SECTIONS values:\n'
                     'goatools.test_data.sections.gjoneska_pfenning \n'
                     'goatools/test_data/sections/gjoneska_pfenning.py \n'
                     'data/gjoneska_pfenning/sections_in.txt\n'))
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
@plac.flg('--prt_study_gos_only',
               help=('Print GO terms only if they are associated with study genes. '
                     'This is useful if printng all GO results '
                     'regardless of their significance (--pval=1.0). '))
@plac.flg('--indent', help="indent GO terms")
@plac.opt('--obo', type=str,
               help="Specifies location and name of the obo file")
@plac.flg('--no_propagate_counts',
               help="Do not propagate counts to parent terms")
# no -r:   args.relationship == False
# -r seen: args.relationship == True
@plac.flg('--relationship',
               help="Propagate counts up all relationships,")
# NO --relationships                -> None
# --relationships part_of regulates -> relationships=['part_of', 'regulates']
# --relationships=part_of           -> relationships=['part_of']
# --relationships=part_of,regulates -> relationships=['part_of', 'regulates']
# --relationships=part_of regulates -> NOT VALID
@plac.flg('--relationships', nargs='*',
               help=('Propagate counts up user-specified relationships, which include: '
                     '{RELS}').format(RELS=' '.join(RELATIONSHIP_LIST)))

@plac.opt('--method', default="bonferroni,sidak,holm,fdr_bh", type=str,
               help=Methods().getmsg_valid_methods())
@plac.opt('--pvalcalc', default="fisher", type=str,
               help=str(FisherFactory()))
@plac.opt('--min_overlap', default=0.7, type=float,
               help="Check that a minimum amount of study genes are in the population")
@plac.opt('--goslim', default='goslim_generic.obo', type=str,
               help="The GO slim file is used when grouping GO terms.")
@plac.opt('--ev_inc', type=str,
               help="Include specified evidence codes and groups separated by commas")
@plac.opt('--ev_exc', type=str,
               help="Exclude specified evidence codes and groups separated by commas")
@plac.flg('--ev_help', help="Print all Evidence codes, with descriptions")
@plac.flg('--ev_help_short', help="Print all Evidence codes")
def run(*taxids):

    # Set the verbosity
    utils.set_verbosity(logger, level=int(0))

    # Load study, population, associations, and GoDag. Run GOEA.
    obj = GoeaCliFnc(GoeaCliArgs().args)

    # Reduce results to significant results (pval<value)
    results_specified = obj.get_results()

    # Print results in a flat list
    obj.prt_results(results_specified)

    pass

