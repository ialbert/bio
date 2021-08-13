import os
import string
import sys
from itertools import *

import plac
from Bio import SeqIO, Seq, SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from biorun import utils

DNA, PEP = "DNA", "PEP"

LOCAL_ALIGN, GLOBAL_ALIGN, SEMIGLOBAL_ALIGN = 1, 2, 3

NUCLEOTIDE, PEPTIDE = "nucleotide", "peptide"

NUCS = set("ATGC" + 'atgc')
PEPS = set("ACDEFGHIKLMNPQRSTVWY")
RNAS = set("AUGC" + 'augc')


def is_nuc(c):
    return c in NUCS


def is_prot(c):
    return c in PEPS


def all_pep(text, limit=1000):
    subs = text[:limit]
    return all(map(is_prot, subs))


def all_nuc(text, limit=1000):
    subs = text[:limit]
    return all(map(is_nuc, subs))


def maybe_sequence(text):
    if all_nuc(text):
        return NUCLEOTIDE
    elif all_pep(text):
        return PEPTIDE
    else:
        return None


def parse(text, idx=0):
    # Perhaps it is already a stream.
    if hasattr(text, "read", ):
        recs = list(SeqIO.parse(text, format='fasta'))
    elif os.path.isfile(text):
        stream = open(text)
        recs = list(SeqIO.parse(stream, format='fasta'))
    elif maybe_sequence(text):
        seq = Seq.Seq(text)
        sid = f"Seq{idx}"
        rec = SeqRecord.SeqRecord(seq=seq, id=sid, name=sid, description='')
        recs = [rec]
    else:
        utils.error(f"Invalid file/sequence: {text}")
        recs = []

    return recs


class Param():
    """Placeholder for parameters"""

    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def print_trace(par, width=80):
    """
    Prints an alignment trace
    """

    for start in range(0, len(par.trace), width):
        end = min(start + width, par.tlen)
        seq1 = par.seqA[start:end]
        trcX = par.trace[start:end]
        seq2 = par.seqB[start:end]

        print(seq1)
        print(trcX)
        print(seq2)
        print()


def format_aln(target, query, aln, par):
    """
    Returns an object with alignment information all set.
    """
    seqA, trace, seqB = format(aln).splitlines()

    # Find non empty indices in the trace
    indices = list(filter(lambda x: x[1] != ' ', enumerate(trace)))
    start, end = indices[0][0], indices[-1][0] + 1

    # Alter the trace to make it look nicer
    trace = trace.replace("-", " ")

    # Populate alignment results
    fmt = Param(**par.__dict__)
    fmt.target = target
    fmt.query = query
    fmt.aln = aln
    fmt.par = par
    fmt.seqA = seqA[start:end]
    fmt.trace = trace[start:end]
    fmt.seqB = seqB[start:end]
    fmt.tlen = len(trace)
    fmt.alen = len(seqA)
    fmt.blen = len(seqB)
    fmt.ident = trace.count('|')
    fmt.mis = trace.count('.')
    fmt.dels = seqB.count("-")
    fmt.ins = seqA.count("-")
    fmt.gap = fmt.dels + fmt.ins
    fmt.pident = fmt.ident / fmt.blen * 100

    return fmt


def print_default(fmt):
    print()
    if fmt.is_dna:
        label = "DNA"
    else:
        label = "PEP"

    print(f"# {label}: {fmt.target.id} ({len(fmt.target):,}) vs {fmt.query.id} ({len(fmt.query):,})")
    print(
        f"# Length={fmt.tlen} Ident={fmt.ident}/{fmt.blen}({fmt.pident:0.1f}%) Mis={fmt.mis} Del={fmt.dels} Ins={fmt.ins} Gap={fmt.gap}")

    if fmt.matrix:
        print(f"# Matrix={fmt.matrix}", end=' ')
    else:
        print(f"# Match={fmt.match} Mismatch={fmt.mismatch}", end=' ')

    print (f"Open={fmt.gap_open} Extend={fmt.gap_extend}")
    print()

    print_trace(fmt)


def align(target, query, par):

    # Query and target sequences.
    t = str(target.seq).upper()
    q = str(query.seq).upper()

    aligner = PairwiseAligner()

    # Select local mode. Global, semiglobal are about scoring.
    if par.mode == LOCAL_ALIGN:
        aligner.mode = 'local'
    else:
        aligner.mode = 'global'

    # Attempts to detect DNA vs peptide sequences.
    if not par.type:
        par.is_dna = all_nuc(t) and all_nuc(q)
    else:
        par.is_dna = par.type == DNA

    aligner.match_score = par.match
    aligner.mismatch_score = par.mismatch

    # Default alignment matrix for peptides.
    if not par.is_dna and not par.matrix:
        par.matrix = 'BLOSUM62'

    # Attempt to load the matrix if specified.
    if par.matrix:
        aligner.substitution_matrix = get_matrix(par.matrix)

    # Internal gap scoring.
    aligner.open_gap_score = -par.gap_open
    aligner.extend_gap_score = -par.gap_extend

    # Global alignment.
    if par.mode == GLOBAL_ALIGN:
        aligner.target_end_open_gap_score = -par.gap_open
        aligner.target_end_extend_gap_score = -par.gap_extend

    # Semiglobal alignment.
    if par.mode == SEMIGLOBAL_ALIGN:
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

    # Performs the alignment
    alns = aligner.align(t, q)

    return alns



    alns = map(builder, alns)

    # Format the aligners
    if par.table:
        print_func = print_tabular
    elif par.mutations:
        print_func = print_mutations
    else:
        print_func = print_pairwise

    for index, aln in enumerate(alns):
        print_func(aln, param=par, index=index)

def get_matrix(matrix, show=False):
    try:

        if os.path.isfile(matrix):
            if show:
                print(open(matrix).read(), end='')
            mat = substitution_matrices.read(matrix)
        else:
            mat = substitution_matrices.load(matrix)
            if show:
                path = os.path.dirname(os.path.realpath(substitution_matrices.__file__))
                path = os.path.join(path, "data", matrix)
                print(open(path).read(), end='')
    except Exception as exc:
        valid = substitution_matrices.load()
        utils.error(f"error loading matrix: {matrix}", stop=False)
        utils.error(f"{exc}", stop=False)
        utils.error(f"valid values: {', '.join(valid)}")


    return mat

@plac.pos("sequence", "sequences")
@plac.opt("match", "match", type=int, abbrev='m')
@plac.opt("mismatch", "mismatch", type=int, abbrev='s')
@plac.opt("gap_open", "gap_open", type=int, abbrev='o')
@plac.opt("gap_extend", "gap_extend", type=int, abbrev='x')
@plac.opt("matrix", "matrix", abbrev='M')
@plac.flg("local_", "local alignment", abbrev='L')
@plac.flg("global_", "local alignment", abbrev='G')
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.opt("type_", "sequence type (nuc, pep)", choices=[DNA, PEP])
def run(gap_open=11, gap_extend=1, matrix='', match=5, mismatch=4, local_=False, global_=False,
        semiglobal=False, type_='',  *sequences):

    # Keeps track of the alignment parameters.
    par = Param()
    par.matrix = None
    par.gap_open = gap_open
    par.gap_extend = gap_extend
    par.mode = SEMIGLOBAL_ALIGN
    par.matrix = matrix
    par.match = match
    par.mismatch = mismatch
    par.type = type_

    # Select alignment mode.
    if local_:
        par.mode = LOCAL_ALIGN
    elif global_:
        par.mode = GLOBAL_ALIGN
    elif semiglobal:
        par.mode = SEMIGLOBAL_ALIGN

    # Input data sources.
    lines = []

    if not sys.stdin.isatty():
        lines.append(sys.stdin)

    # Command line will be second in line.
    lines.extend(sequences)

    counter = cycle(string.ascii_uppercase)

    recs = []
    for idx, text in zip(counter, lines):
        recs.extend(parse(text, idx=idx))

    # If only matrix is specified print it to the screen.
    if matrix and len(recs) == 0:
        get_matrix(matrix, show=True)
        sys.exit()

    # Sequences must be present to be aligned.
    if len(recs) < 2:
        utils.error("Need at least two sequences to align")

    # Keeping people from accidentally running alignments that are too large.
    MAXLEN = 30000
    for rec in recs:
        if len(rec) > MAXLEN:
            utils.error("We recommend that you use a different software.", stop=False)
            utils.error(f"Sequence {rec.id} is too long for this aligner: {len(rec)} > MAXLEN={MAXLEN:,}")

    target = recs[0]

    for query in recs[1:]:

        alns = align(target, query, par=par)

        for aln in alns:
            fmt = format_aln(target=target, query=query, aln=aln, par=par)

            print_default(fmt)

            break


if __name__ == '__main__':
    run()
