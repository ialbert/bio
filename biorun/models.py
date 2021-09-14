from . import utils
from itertools import *


class Sequence:
    """
    Represents a sequence with annotations.
    """

    def __init__(self, title, seq='', type='', ann={}):

        elems = title.strip().split(maxsplit=2)
        if len(elems) == 2:
            self.name, self.desc = elems
        else:
            self.name, self.desc = elems[0], ''
        self.type = type
        self.ann = ann
        self.seq = seq.upper()


class Alignment:
    """
    Represents a pairwise alignment.
    """

    def __init__(self, target: Sequence, query: Sequence, score, par, trace='', **kwds):

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

        if trace:
            # Figuring out start/end padding from a BioPython trace
            pad_func = lambda x: x == ' '
            pad_count = lambda x: sum(1 for _ in takewhile(pad_func, x))

            start = pad_count(trace)
            chop = pad_count(reversed(trace))
            end = len(self.target.seq) - chop

            # Slice only if needed
            if start or chop:
                self.target.seq = self.target.seq[start:end]
                self.query.seq = self.query.seq[start:end]

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

        # Percent identity (BLAST identity)
        # Target sequence representes a gapped alignment.
        # Alternative ways to compute it: https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
        self.pident = self.ident / (self.ident + self.mis + self.dels + self.ins) * 100
        self.alen = len(self.target.seq.strip("-"))

        # Override defaults
        self.par = par
        self.__dict__.update(kwds)


def format_pairwise(aln):
    """
    Format an alignment in pairwise mode
    """
    label = 'DNA' if aln.is_dna else 'PEP'
    out = [

        f"# {label}: {aln.target.name} ({aln.tlen:,}) vs {aln.query.name} ({aln.qlen:,}) score={aln.score}",
        f"# Alignment: pident={aln.pident:0.1f}% len={aln.alen} ident={aln.ident} mis={aln.mis} del={aln.dels} ins={aln.ins}",

    ]

    if aln.par.matrix:
        start = f"# Parameters: matrix={aln.par.matrix}"
    else:
        start = f"# Parameters: match={aln.par.match} penalty={aln.par.mismatch}"

    elem = f"{start} gapopen={aln.par.gap_open} gapextend={aln.par.gap_extend}"

    out.append(elem)

    print("\n".join(out))
