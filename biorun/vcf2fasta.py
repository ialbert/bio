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

def parse_json(fname):
    if fname:
        data = json.load(open(fname))
    else:
        data = []

    remap = {}
    for elem in data:
        key = elem.get('run_accession') or 'xyz'
        alias = elem.get('sample_alias') or 'sample'
        date = elem.get('collection_date') or 'date'
        value = f'{date}-{alias}-{key}'
        remap[key] = value

    return remap

@plac.pos('fname', "input vcf file", type=Path)
@plac.opt('json_', "json file to rename the sequences")
def run(fname, json_=''):
    """
    Creates a FASTA representation from a multisample VCF file where each sample corresponds to a FASTA
    record containing the first ALT if present otherwise the REF.
    """
    vcf = VCF(fname, strict_gt=True)

    samples = vcf.samples

    stream = islice(vcf, None)

    #Load the metdata
    remap = parse_json(json_)

    # Stores the sequences for each sample.
    store = []
    for samp in samples:
        key = remap.get(samp, samp)
        store.append((key, []))

    for v in stream:
        # Go over each sample
        combined = zip(store, v.genotypes)

        for ((key, values), gtypes) in combined:
            # First alternate genotype
            for idx in gtypes[:-1]:
                if idx > 0:
                    values.append(v.ALT[idx - 1])
                    break
            else:
                # Executes if there are no variants.
                values.append(v.REF)

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
