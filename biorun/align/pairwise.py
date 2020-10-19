import warnings, sys, os
import plac
from pprint import pprint
from itertools import islice, count
import textwrap

from biorun.data import fetch
from biorun.const import *
from biorun import models
from biorun import utils

try:
    from Bio import SeqIO
    from Bio import Align
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio import BiopythonExperimentalWarning
except ImportError as exc:
    print(f"*** Error: {exc}", file=sys.stderr)
    print(f"*** Please install biopython: conda install -y biopython==1.76", file=sys.stderr)
    sys.exit(-1)

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio.Align import substitution_matrices


def unpack(aln):
    """
    Unpack the formatted alignment.
    """
    query, pattern, target = format(aln).splitlines()

    return query, pattern, target


def print_aln(aln, matrix, query, target, aligner, width=100):

    nw = 8
    tgt_name = f"{target.name[:nw]:8s}"
    pat_name = " " * nw
    rec_name = f"{query.name[:nw]:8s}"

    query, pattern, target = unpack(aln)

    aln_len = len(pattern)

    ident = pattern.count("|")
    ident_perc = round(100 * ident / aln_len, 1)

    gapn = pattern.count('-')
    gapn_perc = round(100 * gapn / aln_len, 1)

    print(f"# Lenght: {aln_len}")
    print(f"# Identity: {ident}/{aln_len} ({ident_perc}%)")
    print(f"# Gaps: {gapn}/{aln_len} ({gapn_perc}%)")
    print(f"# Score: {aln.score}")
    print(f"#")
    print(f"# Matrix: {matrix} ")
    print(f"# Score: {aln.score}")
    print(f"# Gap open: {aligner.internal_open_gap_score}")
    print(f"# Gap extend: {aligner.internal_extend_gap_score}")
    print("#")

    for start in range(0, aln_len, width):
        end = start + width

        print(tgt_name, target[start:end])
        print(pat_name, pattern[start:end])
        print(rec_name, query[start:end])
        print("")

def biopython_align(query, target, nucl=True, gap_open=None, gap_extend=None, matrix=None, limit=1):
    """
    Perform alignment with BioPython.
    """
    import parasail

    # The pairwise aligner.
    aligner = Align.PairwiseAligner()

    # The default scoring matrices
    if nucl:
        gap_open = gap_open or -10
        gap_extend = gap_extend or -0.5
        matrix = matrix or "NUC.4.4"
    else:
        gap_open = gap_open or -16
        gap_extend = gap_extend or -4
        matrix = matrix or "BLOSUM62"

    # Read the substitution matrix
    try:
        if os.path.isfile(matrix):
            m = substitution_matrices.read(matrix)
        else:
            m = substitution_matrices.load(matrix)
    except Exception as exc:
        print(f"*** Unable to read scoring matrix: {exc}", file=sys.stderr)
        print(f"*** Builtin: {', '.join(substitution_matrices.load())}", file=sys.stderr)
        sys.exit(-1)

    # Gap open
    aligner.open_gap_score = gap_open

    # Gap extend
    aligner.extend_gap_score = gap_extend

    # Assign the matrix.
    aligner.substitution_matrix = m

    # Gaps at the end of the sequences.
    aligner.left_gap_score = 0
    aligner.right_gap_score = 0

    seq_q = str(query.seq).upper()
    seq_t = str(target.seq).upper()


    try:
        results = aligner.align(seq_t, seq_q)
    except Exception as exc:
        print(f"*** Error: {exc}", file=sys.stderr)
        sys.exit(-1)

    # How many alignments to report
    results = islice(results, limit)

    for aln in results:
        print_aln(aln, matrix=matrix, query=query, target=target, aligner=aligner)


def parasail_align(query, target, gap_open=10, gap_extend=0.5, matrix=None, limit=1):
    import parasail

    q = str(query.seq)
    t = str(target.seq)

    result = parasail.sw_trace_striped_16(q, t, 11, 1, parasail.nuc44)

    print(dir(result))

    words = "score similar cigar end_ref"
    for word in words.split():
        if hasattr(result, word):
            value = getattr(result, word)
            print(f"{word}={value}")

    if hasattr(result, "cigar"):
        cigar = result.cigar

        # cigars have seq, len, beg_query, and beg_ref properties
        # the seq property is encoded
        print(cigar.seq)
        # use decode attribute to return a decoded cigar string
        print(cigar.decode)

@plac.opt('start', "start coordinate ", type=int)
@plac.opt('end', "end coordinate", type=int)
@plac.opt('matrix', "scoring matrix")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(start=1, end=None, matrix='', verbose=False, query='', target=''):
    "Prints the effect of an annotation"

    # Move to zero based coordinate system.
    start = utils.shift_start(start)

    queries = fetch.get_data(query)
    targets = fetch.get_data(target)

    param = utils.Param(start=start, end=end)

    query = models.get_origin(queries[0], param)[:200]
    target = models.get_origin(targets[0], param)[:200]

    #print (seq1)
    #print (seq2)

    #biopython_align(query=query, target=target, matrix=matrix)

    parasail_align(query=query, target=target)


def main():
    plac.call(run)


if __name__ == '__main__':
    main()
