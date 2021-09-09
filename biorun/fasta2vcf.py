import plac
from operator import itemgetter
from itertools import *

MATCH, SNP, INS, DEL = 'M', 'X', 'INS', 'DEL'


class Record:
    def __init__(self, name, lines):
        self.name = name.rstrip()
        self.seq = "".join(lines).replace(" ", "").replace("\r", "").upper()


def fast_parser(stream):
    """
    Inspired by Bio.SeqIO.FastaIO.SimpleFastaParser
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in stream:
        if line[0] == ">":
            title = line[1:]
            break
    else:
        return

    lines = []
    for line in stream:
        if line[0] == ">":
            yield Record(name=title, lines=lines)
            lines = []
            title = line[1:]
            continue
        lines.append(line.rstrip())

    fasta = Record(name=title, lines=lines)
    yield fasta


def find_variants(ref, tgt):
    stream = zip(ref.seq, tgt.seq)

    # stream = islice(stream, 50000)

    collect = []
    variants = []
    pos = 0
    lastop = None
    for a, b in stream:

        # Not a valid pairwise alignment
        if a == '-' and b == '-':
            raise Exception(f"Input not a pairwise alignment: {a} vs {b}")

        # Matching bases.
        if a == b:
            pos += 1
            op = MATCH
        elif b == '-':
            # Deletion from query.
            pos += 1
            op = DEL
        elif a == '-':
            op = INS
        elif a != b:
            # Mismatching bases.
            pos += 1
            op = SNP
        else:
            raise Exception(f"Input not a pairwise alignment: {a} vs {b}")

        if lastop != op:
            if collect:
                if lastop != MATCH:
                    variants.append((lastop, collect))
            # Reset the collector
            collect = []

        collect.append((pos, a, b))

        lastop = op

    if collect:
        variants.append((lastop, collect))

    return variants


def format_variants(ref, tgt, variants):
    # This is necessary because consecutive variants may overlap SNP + INSERT for example
    # The later variant wins
    vardict = dict()

    for key, elems in variants:

        if key == SNP:
            for pos, base, alt in elems:
                uid = f"{pos}_{base}_{alt}"
                value = [ref.name, str(pos), uid, base, alt, ".", "PASS", ".", "GT", "1"]
                vardict[pos] = value

        elif key == DEL or key == INS:
            pos = elems[0][0]
            base = ''.join(e[1] for e in elems).strip("-")
            alt = ''.join(e[2] for e in elems).strip("-")

            if pos > 1:
                pos = pos
                base = ref.seq[pos - 1] + base
                alt = tgt.seq[pos - 1] + alt
            uid = f"{pos}_{key}_{len(base + alt)}"
            value = [ref.name, str(pos), uid, base, alt, ".", "PASS", ".", "GT", "1"]
            vardict[pos] = value

    return vardict


def print_variants(ref, tgt, vardict):
    print('##fileformat=VCFv4.2')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FILTER=<ID=PASS,Description="All filters passed">')
    print(f'##contig=<ID={ref.name},length={len(ref.seq.strip("-"))},assembly={ref.name}>')
    print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{tgt.name}")

    for value in vardict.values():
        print("\t".join(value))


@plac.pos("fasta")
def run(*fname):
    fname = '../test/data/mafft.fa'
    stream = open(fname)

    recs = fast_parser(stream)

    ref = next(recs)
    tgt = next(recs)

    variants = find_variants(ref, tgt)

    vardict = format_variants(ref=ref, tgt=tgt, variants=variants)

    print_variants(ref=ref, tgt=tgt, vardict=vardict)


if __name__ == '__main__':
    plac.call(run)
