"""
Handles FASTA related outputs
"""

import biorun.libs.placlib as plac
from biorun.models import jsonrec
from biorun import utils, const, fetch, objects

logger = utils.logger


def get_fasta(item, param):
    """
    Returns a fasta seqrecord for a JSON item.
    """

    # Ignore translation warnings when the sequence in not a multiple of 3.
    if param.translate:
        import warnings
        from Bio import BiopythonWarning
        warnings.simplefilter('ignore', BiopythonWarning)

    if param.inter:
        # Interactive mode returns the origin.
        recs = jsonrec.get_origin(item, param=param)
    elif param.protein:
        recs = jsonrec.get_translation_records(item, param=param)
    elif param.features:
        recs = jsonrec.get_feature_records(item, param=param)
    else:
        recs = jsonrec.get_origin(item, param=param)

    return recs


def print_fasta(recs):
    """
    Prints fasta records
    """
    for rec in recs:
        print(rec.format("fasta"))


def fasta_view(params):
    """
    Converts data to fastya
    """

    for param in params:

        # Stop when data was not found.
        if not param.json:
            utils.error(f"data not found: {param.acc}")

        # Each data may have multiple entries.
        for item in param.json:

            # Get the fasta for each entry.
            recs = get_fasta(item, param=param)

            # Print the fasta records.
            print_fasta(recs)

