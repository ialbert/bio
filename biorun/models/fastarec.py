"""
Generates FASTA outputs from a JSON record.
"""
from biorun import utils
from biorun.models import jsonrec

logger = utils.logger


def get_fasta(item, param):
    """
    Identifies which type of data to return based on parameters.
    """

    # If there is no other filtering, produce the origin.
    origin = not (param.gene or param.type or param.protein or param.translate)

    if origin:
        recs = jsonrec.get_origin(item, param=param)
    elif param.protein:
        recs = jsonrec.get_translation_records(item, param=param)
    else:
        recs = jsonrec.get_feature_records(item, param=param)

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
            utils.error(f"data not found: {param.name}")

        # Each data may have multiple entries.
        for item in param.json:

            # Get the fasta for each entry.
            recs = get_fasta(item, param=param)

            # Print the fasta records.
            print_fasta(recs)
