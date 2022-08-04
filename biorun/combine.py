"""
This module combines kallisto and salmon abundance files scattered across directories.

The output of the program is similar to the count file produced by featureCounts

Usage:
    
     bio combine DIRS
    
Where:

 - DIRS are directories containing

For each id the IDS file the program will attempt to load the file named as 

   DIR/abundances.tsv or DIR/quant.sf files.

Once all files are loaded these are combined into a result that has six columns to start each row, 
then the estimated abundance from each file labeled with a sample name.
"""
import csv
import os
import sys

from biorun.libs import placlib as plac


def parse(dirname):
    """
    Returns a generator over a file
    """

    cands = ["abundance.tsv", "quant.sf"]
    cands = list(map(lambda x: os.path.join(dirname, x), cands))
    exist = list(map(os.path.isfile, cands))

    if not any(exist):
        print("# Not found: %s" % " or ".join(cands))
        sys.exit()

    path = cands[0] if exist[0] else cands[1]

    # Open the stream.
    stream = csv.DictReader(open(path), delimiter="\t")

    # Generate the data
    for row in stream:
        yield row


@plac.pos("names", "input directory names")
def run(*names):
    """
    Processes and combines tsv files.
    """

    names = list(names)

    for name in names:
        if not os.path.isdir(name):
            print("# *** directory not found: %s" % name, file=sys.stderr)
            sys.exit()


    collect = {}

    for dirname in names:

        rows = list(parse(dirname))

        for row in rows:

            try:
                target_id, length, eff_length = row['target_id'], row['length'], row['eff_length']
            except KeyError as exc:
                target_id, length, eff_length = row['Name'], row['Length'], row['EffectiveLength']

            # Fill in the first six columns.
            if target_id not in collect:
                collect[target_id] = [target_id, length, eff_length, "0", "0", "0"]

            try:
                value = row['est_counts']
            except KeyError as exc:
                value = row['NumReads']

            value = str(round(float(value), 1))
            collect[target_id].append(value)

    # Print the collected data
    header = ["target_id", "length", "eff_length", "X1", "X2", "X3"] + names
    header = "\t".join(header)
    print(header)

    # Print the header

    for key, values in collect.items():
        line = "\t".join(values)
        print(line)


if __name__ == '__main__':
    plac.call(run)
