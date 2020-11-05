"""
Generates GFF outputs from a JSON record.
"""



def gff_view(params):
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
