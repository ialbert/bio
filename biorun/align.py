import os
from itertools import *

import plac, string
from Bio import SeqIO, Seq, SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from biorun import utils

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
    if os.path.isfile(text):
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
    pass


def print_trace(seqA, seqB, trace, width=80):
    """
    Prints an alignment trace
    """
    alen = len(trace)

    for start in range(0, len(trace), width):
        end = min(start + width, alen)
        seq1 = seqA[start:end]
        trcX = trace[start:end]
        seq2 = seqB[start:end]

        print(seq1)
        print(trcX)
        print(seq2)
        print()


def trace_counter(trace, mchr='.', ichr='|', dchr='-'):
    count_mis = trace.count(mchr)
    count_ident = trace.count(ichr)
    return count_ident, count_mis


def print_aln(target, query, aln):
    seqA, trace, seqB = format(aln).splitlines()

    # Find non empty indices in the trace
    indices = list(filter(lambda x: x[1] != ' ', enumerate(trace)))
    start, end = indices[0][0], indices[-1][0] + 1

    seqA = seqA[start:end]
    trace = trace[start:end]
    seqB = seqB[start:end]

    # Alter the trace
    trace = trace.replace("-", " ")

    alen = len(trace)
    blen = len(seqB)
    count_ident, count_mis = trace_counter(trace)
    count_delete = seqB.count("-")
    count_insert = seqA.count("-")
    count_gap = count_delete + count_insert
    perc_idn = count_ident / len(seqB) * 100


    print()
    print(f"# {target.id} ({len(target):,}) vs {query.id} ({len(query):,})")
    print(f"# Length={alen} Ident={count_ident}/{blen}({perc_idn:0.1f}%) Mis={count_mis} Del={count_delete} Ins={count_insert} Gap={count_gap}")
    print()
    print_trace(seqA, seqB, trace)


def align(target, query, param):
    # Query and target sequences.
    t = str(target.seq).upper()
    q = str(query.seq).upper()

    aligner = PairwiseAligner()

    # Select local mode. Global, semiglobal are about scoring.
    if param.mode == LOCAL_ALIGN:
        aligner.mode = 'local'

    # Attempts to detect DNA vs peptide sequences.
    param.is_dna = all_nuc(t) and all_nuc(q)

    # Default substituion matrix.
    if not param.matrix:
        param.matrix = 'NUC.4.4' if param.is_dna else 'BLOSUM62'

    # Apply substitution matrix.
    aligner.substitution_matrix = substitution_matrices.load(param.matrix)

    # Internal gap scoring.
    aligner.open_gap_score = -param.gap_open
    aligner.extend_gap_score = -param.gap_extend

    # Global alignment.
    if param.mode == GLOBAL_ALIGN:
        aligner.target_end_open_gap_score = -param.gap_open
        aligner.target_end_extend_gap_score = -param.gap_extend

    # Semiglobal alignment.
    if param.mode == SEMIGLOBAL_ALIGN:
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

    # Performs the alignment
    alns = aligner.align(t, q)

    return alns

    # Reformat alignments as a more detailed class.
    def builder(aln):
        rec = Alignment(qseq=qseq, tseq=tseq, aln=aln, param=param)
        return rec

    alns = map(builder, alns)

    # Format the aligners
    if param.table:
        print_func = print_tabular
    elif param.mutations:
        print_func = print_mutations
    else:
        print_func = print_pairwise

    for index, aln in enumerate(alns):
        print_func(aln, param=param, index=index)


@plac.pos("sequence", "sequences")
@plac.opt("gap_open", "gap_open", type=int, abbrev='o')
@plac.opt("gap_extend", "gap_extend", type=int, abbrev='e')
@plac.flg("local_", "local alignment", abbrev='L')
@plac.flg("global_", "local alignment", abbrev='G')
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.opt("type_", "sequence type (nuc, pep)", choices=["nuc", "pep"])
def run(type_="nuc", gap_open=11, gap_extend=1, local_=False, global_=False, semiglobal=False, *sequences):
    param = Param()
    param.matrix = None
    param.gap_open = gap_open
    param.gap_extend = gap_extend
    param.mode = SEMIGLOBAL_ALIGN

    if local_:
        param.mode = LOCAL_ALIGN
    elif global_:
        param.mode = GLOBAL_ALIGN
    elif semiglobal:
        param.mode = SEMIGLOBAL_ALIGN

    lines = []

    # if not sys.stdin.isatty():
    #    lines.extend(sys.stdin)

    lines.extend(sequences)

    counter = cycle(string.ascii_uppercase)

    recs = []
    for idx, text in zip(counter, lines):
        recs.extend(parse(text, idx=idx))

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
        alns = align(target, query, param=param)

        for aln in alns:
            print_aln(target, query, aln)
            break


if __name__ == '__main__':
    run()
