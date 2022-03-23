import json
import sys
from pathlib import Path

import biorun.libs.placlib as plac

try:
    from cyvcf2 import VCF
except ImportError as exc:
    print(f"# Error: {exc}")
    print(f"# Use: pip install cyvcf2")
    sys.exit(1)

from itertools import *

LIMIT = None


def parse_json(fname, patt="collection_date:sample_alias", key_name="run_accession"):

    # Individual patterns to match
    patts = patt.split(":")

    # Load the json file.
    if fname:
        data = json.load(open(fname))
    else:
        data = []

    # Remapping dictionary.
    remap = {}
    counter = count(1)
    for idx, elem in zip(counter, data):

        # Insert the index into the dictionary.
        elem['idx'] = str(idx)

        # This field must be present to rename.
        try:
            key = elem[key_name]
        except KeyError as exc:
            sys.stderr.write(f"# Metadata error in: {fname} \n")
            sys.stderr.write(f"# Key not found: {key_name}\n")
            sys.exit(1)

        # Build the rename pattern
        fields = []
        for p in patts:
            fields.append(elem.get(p, p))

        # The rename value.
        value = "_".join(fields)

        # Last pattern is the key.
        remap[key] = value

    return remap


@plac.pos('fname', "input vcf file", type=Path)
@plac.opt('json_', "json file to rename the sequences")
@plac.opt('chrom', "chromomosome to process (default: first sequence)")
@plac.opt('limit', "how many variants to process (0=no limit)", type=int)
@plac.opt('key', "renaming key (advanced use only)")
@plac.opt('patt', "renaming pattern (advanced use only)")
@plac.flg('ref', "generate a reference as well")
def run(fname, json_="", chrom="", limit=0, key="", patt="", ref=False):
    """
    Creates a FASTA representation from a multisample VCF file where each sample corresponds to a FASTA
    record containing the first ALT if present otherwise the REF.
    """

    # Renaming key.
    key_name= key or "run_accession"

    # Renaming pattern
    patt = patt or "collection_date:sample_alias:run_accession"

    # Parse the VCF file.
    vcf = VCF(fname, strict_gt=True)

    # Select the target chromosome.
    chrom = chrom if chrom else vcf.seqnames[0]

    # Chromosome filtering function.
    def chrom_filter(v):
        return v.CHROM == chrom

    # Load the metdata
    remap = parse_json(json_, patt=patt, key_name=key_name)

    # Iterate over the VCF records.

    # Stores the reference sequence
    refseq = []

    # Stores the sequences for each sample.
    store = []

    # Initialize the sequence list.
    for samp in vcf.samples:
        key = remap.get(samp, samp)
        store.append((key, []))

    # Set the limits.
    limit = None if limit < 1 else limit

    # Used during debugging mostly.
    stream = islice(vcf, limit)

    # Apply the chromosome filter.
    stream = filter(chrom_filter, stream)

    # Process each variant
    for v in stream:

        # Parallel iteration of sequences and genomtypes.
        combined = zip(store, v.genotypes)

        # All valid bases at the current index
        bases = [v.REF] + v.ALT

        # Collect the reference sequence.
        refseq.append(v.REF)

        # Select the first ALT in each genoytpe.
        for ((key, values), gtypes) in combined:
            idx = gtypes[0]
            idx = 0 if idx < 0 else idx
            base = bases[idx]
            values.append(base)

    # Print the reference sequence.
    if ref:
        seq = "".join(refseq)
        print(f">reference\n{seq}")

    # Output each individual sequence.
    for (key, values) in store:
        seq = "".join(values)
        print(f">{key}\n{seq}")


def main():
    """
    Entry point for the script.
    """
    if len(sys.argv) < 2:
        sys.argv.append("-h")
    plac.call(run)


if __name__ == '__main__':
    main()
