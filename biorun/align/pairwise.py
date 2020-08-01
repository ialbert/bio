import warnings
import sys
import plac
from pprint import pprint
from itertools import islice, count
import textwrap

import os

try:
    from Bio import SeqIO
    from Bio import Align
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio import BiopythonExperimentalWarning
except ImportError as exc:
    print(f"*** Error: {exc}", file=sys.stderr)
    print(f"*** Please install biopython: conda install -y biopython", file=sys.stderr)
    sys.exit(-1)

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio.Align import substitution_matrices


def guess_type(name):
    if name.endswith("gb") or name.endswith("genbank"):
        return "genbank"
    return "fasta"


def get_sequence(fname, name, nucl=False):

    if os.path.isfile(fname):
        ftype = guess_type(fname)
        recs = SeqIO.parse(fname, ftype)
    else:
        alph = IUPAC.ambiguous_dna if nucl else IUPAC.protein
        seq = Seq(fname, alphabet=alph)
        rec = SeqRecord(seq, id=name, name=name)
        recs = [rec]
    return recs


def unpack(aln):
    """
    Unpack the formatted alignment.
    """
    query, pattern, target = format(aln).splitlines()

    return query, pattern, target


def print_aln(aln, matrix, tgt, rec, aligner, width=100):

    nw = 8
    tgt_name = f"{tgt.name[:nw]:8s}"
    pat_name = " " * nw
    rec_name = f"{rec.name[:nw]:8s}"

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



def align(query, target, nucl=False, gap_open=0, gap_extend=0, matrix=None, limit=1, report=1, params=None):
    # The query sequences.
    queries = get_sequence(query, name="a")

    # The alignment targets.
    targets = get_sequence(target, name="b")

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


    for tgt in targets:
        tgt_seq = str(tgt.seq).upper()
        # print (tgt)
        for rec in queries:

            rec_seq = str(rec.seq).upper()

            try:
                results = aligner.align(tgt_seq, rec_seq)
            except KeyError as exc:
                print(f"*** Error: {exc}", file=sys.stderr)
                sys.exit(-1)

            # How many alignments to report
            results = islice(results, limit)

            for aln in results:
                print_aln(aln, matrix=matrix, tgt=tgt, rec=rec, aligner=aligner)


@plac.annotations(
    target=("target", "positional", None, str, None, None),
    query=("query", "positional", None, str, None, None),
    nucl=("sequence", "flag",),
    matrix=("scoring matrix", "option", "m"),

)
def run(target, query, nucl, matrix):
    "Prints the effect of an annotation"


    align(query=query, target=target, nucl=nucl, matrix=matrix)


def main():
    plac.call(run)


if __name__ == '__main__':
    main()
