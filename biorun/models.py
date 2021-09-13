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
            self.name, self.desc = elems, ''
        self.type = type
        self.ann = ann
        self.seq = seq.upper()


class Alignment:
    """
    Represents a pairwise alignment.
    """

    def __init__(self, target: Sequence, query: Sequence, trace='', **kwds):

        # Setting overrideable defaults
        self.is_dna = True

        # Original query and target.
        self.query = query
        self.target = target

        # Alignment parameters
        self.ident = self.ins = self.dels = self.mism = 0

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
                self.mism += 1

        # Override defaults
        self.__dict__.update(kwds)


def pairwise_format(aln):

    label = 'DNA' if aln.is_dna else 'PEP'
    out = [

        f"# {label}: {aln.target.name} ({len(aln.target.seq):,}) vs {aln.query.name} ({len(aln.query.seq):,}) score = {aln.score}",

        f"# Alignment: pident={aln.pident:0.1f}% len={aln.tlen} ident={aln.ident} mis={aln.mis} del={aln.dels} ins={aln.ins}",

    ]