from . import utils
from itertools import *
from Bio.Seq import Seq


class Sequence:
    """
    Represents a sequence with a name.
    """

    def __init__(self, title='', seq=''):

        elems = title.strip().split(maxsplit=2)

        self.name = self.desc = self.seq = ''

        if len(elems) == 2:
            self.name, self.desc = elems
        elif elems:
            self.name = elems[0]

        self.seq = seq.upper()

    def __str__(self):
        return f">{self.name} {self.seq[:10]}..."


class Alignment:
    """
    Represents a pairwise alignment.
    """

    def __init__(self, target: Sequence, query: Sequence, score, par, **kwds):

        # Make the trace nicer.

        # Setting overrideable defaults
        self.is_dna = True

        # Original query and target.
        self.query = query
        self.target = target

        # Alignment parameters
        self.ident = self.ins = self.dels = self.mis = 0
        self.score = score
        self.tlen = par.tlen
        self.qlen = par.qlen

        # Compute alignment information.
        for idx, a, b in zip(count(), self.query.seq, self.target.seq):
            if a == b:
                self.ident += 1
            elif a == '-':
                self.ins += 1
            elif b == '-':
                self.dels += 1
            elif a != b:
                self.mis += 1

        # trace = trace.replace("-"," ")

        # Percent identity (BLAST identity)
        # Target sequence representes a gapped alignment.
        # Alternative ways to compute it: https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
        self.pident = self.ident / (self.ident + self.mis + self.dels + self.ins) * 100
        self.alen = len(self.target.seq.strip("-"))

        # Override defaults
        self.par = par
        self.__dict__.update(kwds)


class Param:
    """
    Represents alignment parameters
    """

    def __init__(self, **kwds):
        self.is_dna = False
        self.matrix = None
        self.query = Sequence()
        self.target = Sequence()
        self.gap_open = 11
        self.gap_extend = 1
        self.match = 5
        self.mismatch = 4
        self.mode = ''
        self.type = ''
        # Override with defaults
        self.__dict__.update(kwds)

        self.tlen = len(self.target.seq)
        self.qlen = len(self.query.seq)


def format_pairwise(alns, width=81):
    """
    Formats an alignment in pairwise mode
    """
    for aln in alns:
        label = 'DNA' if aln.is_dna else 'PEP'
        out = [
            "",
            f"# {label}: {aln.target.name}({aln.tlen:,}) vs {aln.query.name}({aln.qlen:,}) score={aln.score}",
            f"# Alignment: pident={aln.pident:0.1f}% len={aln.alen} ident={aln.ident} mis={aln.mis} del={aln.dels} ins={aln.ins}",

        ]

        if aln.par.matrix:
            start = f"# Parameters: matrix={aln.par.matrix}"
        else:
            start = f"# Parameters: match={aln.par.match} penalty={aln.par.mismatch}"

        elem = f"{start} gapopen={aln.par.gap_open} gapextend={aln.par.gap_extend}"

        out.extend([elem, ""])

        # Generate the traces
        for start in range(0, len(aln.target.seq), width):
            end = min(start + width, aln.tlen)

            seq1 = aln.target.seq[start:end]
            seq2 = aln.query.seq[start:end]
            trace = ""
            for a, b in zip(seq1, seq2):
                if a == b:
                    trace += '|'
                elif a == '-' or b == '-':
                    trace += '-'
                else:
                    trace += '.'

            out.append(seq1)
            out.append(trace)
            out.append(seq2)

        out.append("")

        print("\n".join(out))


MATCH, SNP, INS, DEL = 'M', 'SNP', 'INS', 'DEL'


def find_variants(ref, tgt):
    """
    Takes two aligned sequences and tabulates what variants can be found at a given position.
    This is the trickiest to get right ... mostly works but may still have bugs.

    Returns a dictionary keyed by position where each position describes the  variant as a tuple.
    """

    stream = zip(count(), ref.seq, tgt.seq)

    # stream = islice(stream, 50000)

    collect = []
    variants = []
    lastop = None
    pos = 0
    for idx, a, b in stream:

        # Multiple sequence alignment
        if a == '-' and b == '-':
            continue

        if a == b:
            # Matches
            pos += 1
            op = MATCH

        elif b == '-':
            # Deletion from query.
            pos += 1
            op = DEL

        elif a == '-':
            # Insertion into query
            op = INS

        elif a != b:
            # Mismatching bases.
            pos += 1
            op = SNP

        else:
            raise Exception(f"Should never hit this (sanity check): {a} vs {b}")

        # Collect variants when the operator changes.
        if lastop != op:
            if collect:
                if lastop != MATCH:
                    variants.append((lastop, collect))
            # Reset the collector
            collect = []

        # Collect the positions that have been visited
        collect.append((idx, pos, a, b))

        # Keep track of the last previous operation
        lastop = op

    # Collect last element
    if collect:
        variants.append((lastop, collect))

    # This is necessary because consecutive variants may overlap SNP + INSERT for example
    vcfdict = dict()

    for key, elems in variants:

        if key == MATCH:
            continue

        # Lenght of variant
        size = len(elems)

        idx, pos = elems[0][0], elems[0][1]

        # Query starts with insertion
        pos = 1 if pos < 1 else pos

        base = ref.seq[idx:idx + size].strip('-')
        alt = tgt.seq[idx:idx + size].strip('-')

        info = f"TYPE={key}"

        if key == SNP:
            # Mismatches printed consecutively
            for idx, pos, base, alt in elems:
                name = f"{pos}{base}/{alt}"
                value = [ref.name, str(pos), name, base, alt, ".", "PASS", info, "GT", "1"]
                vcfdict[pos] = value


        elif key == INS or key == DEL:
            # Handle insertions and deletions.

            # Push back on POS if it is not 1.
            if idx > 0:
                # For insertions the pos does not advance
                if key == DEL:
                    pos = pos - 1
                base = ref.seq[idx - 1] + base
                alt = tgt.seq[idx - 1] + alt
            else:
                # When there is no preceding base the last base must be used
                lastidx = elems[-1][0] + 1
                base = base + ref.seq[lastidx]
                alt = alt + tgt.seq[lastidx]

            alt = alt or '.'
            base = base or '.'

            if key == DEL:
                name = f"{pos}del{len(base[1:])}"
            else:
                name = f"{pos}ins{len(alt[1:])}"

            value = [ref.name, str(pos), name, base, alt, ".", "PASS", info, "GT", "1"]
            vcfdict[pos] = value

    return vcfdict


def format_vcf(alns):
    for aln in alns:
        vcfdict = find_variants(aln.target, aln.query)

        ref, query = aln.target, aln.query

        print('##fileformat=VCFv4.2')
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        print('##FILTER=<ID=PASS,Description="All filters passed">')
        print('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">')
        print(f'##contig=<ID={ref.name},length={len(ref.seq.strip("-"))},assembly={ref.name}>')
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{query.name}")

        for value in vcfdict.values():
            print("\t".join(value))


def format_table(alns, sep="\t"):
    header = "target query score len pident match mism ins del"
    print("\t".join(header.split()))

    for aln in alns:
        data = [
            f"{aln.target.name}", f"{aln.query.name}",
            f"{aln.score}", f"{aln.pident:0.1f}", f"{aln.tlen}",
            f"{aln.ident}", f"{aln.mis}", f"{aln.dels}", f"{aln.ins}"
        ]
        line = sep.join(data)
        print(line)

def format_variants(alns):
    header = "target query pos len type ref alt"
    print("\t".join(header.split()))
    for aln in alns:
        vcfdict = find_variants(aln.target, aln.query)
        for value in vcfdict.values():
            name = value[2]
            pos = value[1]
            base = value[3]
            alt = value[4]
            size = '0'
            info = value[7]
            if "SNP" in info:
                info = 'snp'
                size = 1
            elif 'DEL' in info:
                info = 'del'
                base = base[1:] if pos != '1' else base[:1]
                alt = '-' * len(base)
                size = len(base)

            elif 'INS' in info:
                info = 'ins'
                alt = alt[1:] if pos != '1' else alt[:1]
                base = '-' * len(alt)
                size = len(alt)

            data = [aln.target.name, aln.query.name, pos, str(size), info, base, alt, ]
            print("\t".join(data))