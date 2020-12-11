
from goatools.cli.find_enrichment import GoeaCliArgs
from goatools.cli.find_enrichment import GoeaCliFnc

from biorun import utils
from biorun.libs import placlib as plac

logger = utils.logger


@plac.flg('taxids', "verbose mode, prints more messages")
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

