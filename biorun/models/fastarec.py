"""
Handles FASTA related outputs
"""

import sys, json
import plac
from biorun.models import jsonrec
from biorun import utils, const

logger = utils.logger


def get_records(item, features=False, proteins=False, translate=False):
    """
    Selects the source of the records.
    """

    # Produce the origin by default
    origin = not(features or proteins or translate)

    if origin:
        recs = jsonrec.origin_records(item)
    elif proteins:
        recs = jsonrec.protein_records(item)
    elif features:
        recs = jsonrec.feature_records(item)
    else:
        recs = []

    # Apply the translation here.
    if translate:
        pass
        #def func(f):
        #    expected = first(f, "translation")
        #    # Stop codon is present in the CDS but not in the translation.
        #    observed = str(dna)[:-1]
        #
        #    # Checking for non-standard translations.
        #    if expected and expected != observed:
        #        logger.info(f"translation mismatch for: {rec.id}")

    return recs

@plac.flg("features", "convert the features")
@plac.flg("proteins", "extract embedded protein sequences")
@plac.flg("translate", "translate DNA sequences")
def run(features=False, proteins=False, translate=False):
    """
    Converts data to fastya
    """

    data = json.loads(sys.stdin.read())

    for item in data:

        recs = get_records(item, features=features, proteins=proteins, translate=translate)

        for rec in recs:
            print(rec.format("fasta"), end='')


