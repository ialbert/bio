"""
This package attempts to simplify the BioPython SeqRecords into a simpler, flatter structure that
can be more readily worked with.
"""
import sys, os
from pprint import pprint
from intervaltree import Interval, IntervalTree
from biorun import utils
from functools import lru_cache
from itertools import *

try:
    from Bio import SeqIO
except ImportError as exc:
    print(f"*** {exc}", file=sys.stderr)
    print(f"*** This software requires biopython.", file=sys.stderr)
    print(f"*** Try: conda install biopython", file=sys.stderr)
    sys.exit(1)

logger = utils.logger


def parse_genbank(stream, fmt=utils.GENBANK):
    """
    Parses GenBank into sequence records.
    """
    recs = SeqIO.parse(stream, format=fmt)
    return recs


if __name__ == "__main__":
    import doctest

    doctest.testmod()
