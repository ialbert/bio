import os
import string
import sys
from itertools import *

import plac

try:
    from Bio import SeqIO, Seq, SeqRecord
    from Bio.Align import PairwiseAligner
    from Bio.Align import substitution_matrices
except ImportError as exc:
    print(f"# Error: {exc}", file=sys.stderr)
    print(f"# This program requires biopython", file=sys.stderr)
    print(f"# Install: conda install -y biopython>=1.79", file=sys.stderr)
    sys.exit(-1)

from biorun import utils

DNA, PEP = "DNA", "PEP"

TABLE_FMT, VARIANT_FMT = "table", "variant"

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


def parse(text, counter):
    # Perhaps it is already a stream.
    if hasattr(text, "read", ):
        recs = list(SeqIO.parse(text, format='fasta'))
    elif os.path.isfile(text):
        stream = open(text)
        recs = list(SeqIO.parse(stream, format='fasta'))
    elif maybe_sequence(text):
        idx = next(counter)
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


def print_trace(par, width=81):
    """
    Prints an alignment trace
    """

    for start in range(0, len(par.trace), width):
        end = min(start + width, par.tlen)
        seq1 = par.seqA[start:end]
        trcX = par.trace[start:end]
        seq2 = par.seqB[start:end]

        print(seq1)
        print(trcX, end=' ')
        print(f"{end:,}")
        print(seq2)
        print()


def format_alignment(target, query, aln, par):
    """
    Returns an object with alignment information all set.
    """
    seqA, trace, seqB = format(aln).splitlines()

    # Find non empty indices in the trace
    indices = list(filter(lambda x: x[1] != ' ', enumerate(trace)))
    start, end = indices[0][0], indices[-1][0] + 1

    # Make the trace nicer
    trace = trace.replace("-", " ")

    # Populate alignment results
    fmt = Param(**par.__dict__)
    fmt.target = target
    fmt.query = query
    fmt.aln = aln
    fmt.score = aln.score
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
    fmt.pident = fmt.ident / fmt.tlen * 100

    return fmt


def print_default(fmt):
    print()
    if fmt.is_dna:
        label = "DNA"
    else:
        label = "PEP"

    print(f"# {label}: {fmt.target.id} ({len(fmt.target):,}) vs {fmt.query.id} ({len(fmt.query):,}) score={fmt.score}")
    print(
        f"# Alignment: pident={fmt.pident:0.1f}% len={fmt.tlen} ident={fmt.ident} mis={fmt.mis} del={fmt.dels} ins={fmt.ins}")

    if fmt.matrix:
        print(f"# Parameters: matrix={fmt.matrix}", end=' ')
    else:
        print(f"# Parameters: match={fmt.match} penalty={fmt.mismatch}", end=' ')

    print(f"gapopen={fmt.gap_open} gapextend={fmt.gap_extend}")
    print()

    print_trace(fmt)

def table_fmt(fmt, sep="\t"):

    data = [
        f"{fmt.target.id}", f"{fmt.query.id}",
        f"{fmt.score}", f"{fmt.pident:0.1f}", f"{fmt.tlen}",
        f"{fmt.ident}", f"{fmt.mis}", f"{fmt.dels}", f"{fmt.ins}"
    ]
    line = sep.join(data)
    print(line)


def variant_fmt(fmt):
    counter = count(1)

    def stream():
        return zip(counter, fmt.seqA, fmt.trace, fmt.seqB)



    def display(data1, data2, oper):
        if not (data1 or data2):
            return
        pos = data1[0][0]
        seq1 = ''.join(x[1] for x in data1)
        seq2 = ''.join(x[1] for x in data2)
        if oper != 'match':
            print(f"{pos}\t{oper}\t{len(seq1)}\t{seq1}\t{seq2}")


    # List the variants one per line
    data1, data2, last = [], [], ''
    for i, a, t, b in stream():

        if t == '.':
            oper = 'mis'
        elif t == '|':
            oper = 'match'
        elif b == '-':
            oper = 'del'
        elif a == '-':
            oper = 'ins'
        else:
            raise Exception()

        if oper != last:
            #print (oper, data1)
            display(data1, data2, last)
            data1, data2, last = [], [], ''

        data1.append((i, a))
        data2.append((i, b))

        last = oper

    # Trailing variants.
    display(data1, data2, last)


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
    aligner.mismatch_score = -par.mismatch

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
    elif par.mode == SEMIGLOBAL_ALIGN:
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

    # Performs the alignment
    alns = aligner.align(t, q)

    return alns


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

def print_header(text):
    if text:
        elems = text.split()
        print ("\t".join(elems))


@plac.pos("sequence", "sequences")
@plac.opt("match", "match", type=int, abbrev='m')
@plac.opt("mismatch", "mismatch", type=int, abbrev='s')
@plac.opt("open_", "gap open penalty", type=int, abbrev='o')
@plac.opt("extend", "gap extend penalty", type=int, abbrev='x')
@plac.opt("matrix", "matrix", abbrev='M')
@plac.flg("variant", "output variants only", abbrev='V')
@plac.flg("table", "format output as a table", abbrev="T")
@plac.flg("local_", "local alignment", abbrev='L')
@plac.flg("global_", "local alignment", abbrev='G')
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.opt("type_", "sequence type (nuc, pep)", choices=[DNA, PEP])
def run(open_=6, extend=1, matrix='', match=1, mismatch=2, local_=False, global_=False,
        semiglobal=False, type_='', variant=False, table=False, *sequences):
    # Keeps track of the alignment parameters.
    par = Param()
    par.matrix = None
    par.gap_open = abs(open_)
    par.gap_extend = abs(extend)
    par.mode = SEMIGLOBAL_ALIGN
    par.matrix = matrix
    par.match = abs(match)
    par.mismatch = abs(mismatch)
    par.type = type_
    par.format = ''
    par.showall = False

    if variant:
        par.format = VARIANT_FMT
    elif table:
        par.format = TABLE_FMT

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
    for text in lines:
        recs.extend(parse(text, counter=counter))


    # If only matrix is specified print it to the screen.
    if matrix and len(recs) == 0:
        get_matrix(matrix, show=True)
        sys.exit()

    # Sequences must be present to be aligned.
    if len(recs) < 2:
        utils.error("Need at least two sequences to align")

    # Keeping people from accidentally running alignments that are too large.
    MAXLEN = 50000
    for rec in recs:
        if len(rec) > MAXLEN:
            utils.error("We recommend that you use a different software.", stop=False)
            utils.error(f"Sequence {rec.id} is too long for this aligner: {len(rec)} > MAXLEN={MAXLEN:,}")

    target = recs[0]


    # Select result formatter
    if par.format == TABLE_FMT:
        header = "target query score len pident match mism ins del"
        formatter = table_fmt
    elif par.format == VARIANT_FMT:
        header = "pos type len target query"
        formatter = variant_fmt
    else:
        header = ''
        formatter = print_default

    print_header(header)

    for query in recs[1:]:

        alns = align(target, query, par=par)

        for aln in alns:

            # Format the alignment into a class that carries a multitude of parameters.
            res = format_alignment(target=target, query=query, aln=aln, par=par)

            # Apply the formatter.
            formatter(res)

            if not par.showall:
                break


if __name__ == '__main__':
    run()
