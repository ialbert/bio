import os
import string
import sys
from itertools import *
from . import models

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

TABLE_FMT, VCF_FMT, VARIANTS_FMT, PAIRWISE_FMT = "table", "vcf", 'variants', 'pariwise'

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


def guess_type(text):
    if all_nuc(text):
        return NUCLEOTIDE
    elif all_pep(text):
        return PEPTIDE
    else:
        return None


def parse_input(obj, counter):
    if hasattr(obj, "read", ):
        # Object may already be a stream
        recs = list(utils.fasta_parser(obj))
    elif os.path.isfile(obj):
        # Perhaps it is a filename
        stream = open(obj)
        recs = list(utils.fasta_parser(stream))
    elif guess_type(obj):
        # Perhaps the sequence comes from command line
        idx = next(counter)
        name = f"Seq{idx}"
        recs = [utils.Fasta(name=name, lines=[obj])]
    else:
        utils.error(f"Invalid file/sequence: {obj}")
        recs = []

    return recs


class Param():
    """Placeholder for parameters"""

    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def format_alignment(target, query, aln, par):
    """
    Returns an object with alignment information all set.
    """
    seqA, trace, seqB = format(aln).splitlines()

    t = models.Sequence(title=target.name, seq=seqA)
    q = models.Sequence(title=query.name, seq=seqB)

    # Work in progress

    par.tlen = len(target.seq)
    par.qlen = len(query.seq)

    al = models.Alignment(target=t, query=q, score=aln.score, trace=trace, par=par)

    models.format_pairwise(al)

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
    fmt.seqA = utils.Fasta(name=target.name, seq=seqA[start:end])
    fmt.trace = trace[start:end]
    fmt.seqB = utils.Fasta(name=query.name, seq=seqB[start:end])

    fmt.tlen = len(fmt.trace)
    fmt.alen = len(seqA)
    fmt.blen = len(seqB)
    fmt.ident = trace.count('|')
    fmt.mis = trace.count('.')
    fmt.dels = seqB.count("-")
    fmt.ins = seqA.count("-")
    fmt.start = start
    fmt.end = end
    fmt.gap = fmt.dels + fmt.ins
    fmt.pident = fmt.ident / fmt.tlen * 100

    return fmt


def pairwise_fmt(fmt, width=81):
    print()
    if fmt.is_dna:
        label = "DNA"
    else:
        label = "PEP"

    print(
        f"# {label}: {fmt.target.name} ({len(fmt.target.seq):,}) vs {fmt.query.name} ({len(fmt.query.seq):,}) score = {fmt.score}")
    print(
        f"# Alignment: pident={fmt.pident:0.1f}% len={fmt.tlen} ident={fmt.ident} mis={fmt.mis} del={fmt.dels} ins={fmt.ins}")

    if fmt.matrix:
        print(f"# Parameters: matrix={fmt.matrix}", end=' ')
    else:
        print(f"# Parameters: match={fmt.match} penalty={fmt.mismatch}", end=' ')

    print(f"gapopen={fmt.gap_open} gapextend={fmt.gap_extend}")
    print()

    # Generate the traces
    for start in range(0, len(fmt.trace), width):
        end = min(start + width, fmt.tlen)

        seq1 = fmt.seqA.seq[start:end]
        trcX = fmt.trace[start:end]
        seq2 = fmt.seqB.seq[start:end]

        print(seq1)
        print(trcX)
        print(seq2)
        print()





def safe_abs(value):
    try:
        return abs(float(value))
    except Exception as exc:
        utils.error(f"{exc} {value}")


def align(target, query, par: models.Param):
    # Query and target sequences.
    t = str(target.seq)
    q = str(query.seq)

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

    # Set alignment defaults
    if par.gap_open == '':
        par.gap_open = 6 if par.is_dna else 11

    if par.gap_extend == '':
        par.gap_extend = 1

    if par.match == '':
        par.match = 1

    if par.mismatch == '':
        par.mismatch = 2

    aligner.match_score = safe_abs(par.match)
    aligner.mismatch_score = -safe_abs(par.mismatch)

    # Default alignment matrix for peptides.
    if not par.is_dna and not par.matrix:
        par.matrix = 'BLOSUM62'

    # Attempt to load the matrix if specified.
    if par.matrix:
        aligner.substitution_matrix = get_matrix(par.matrix)

    # Internal gap scoring.
    aligner.open_gap_score = -safe_abs(par.gap_open)
    aligner.extend_gap_score = -safe_abs(par.gap_extend)

    # Global alignment.
    if par.mode == GLOBAL_ALIGN:
        aligner.target_end_open_gap_score = aligner.open_gap_score
        aligner.target_end_extend_gap_score = aligner.extend_gap_score

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

#
# NUC=(1,2,6,1)
# PEP=(BLOSUM62, 11,1)
#
@plac.pos("sequence", "sequences")
@plac.opt("match", "match (default: 1 for DNA, BLOSUM62 for PEP)", abbrev='m')
@plac.opt("mismatch", "mismatch (default: 1 for DNA, BLOSUM62 for PEP)", abbrev='s')
@plac.opt("open_", "gap open penalty (default: 6 for DNA, 11 for PEP)", abbrev='o')
@plac.opt("extend", "gap extend penalty (default: 1 for DNA, 1 for PEP)", abbrev='x')
@plac.opt("matrix", "matrix (BLOSUM62 for peptides)", abbrev='M')
@plac.flg("vcf", "output vcf file", abbrev='V')
@plac.flg("table", "output in tabular format", abbrev="T")
@plac.flg("vars", "output variant columns", abbrev="A")
@plac.flg("local_", "local alignment", abbrev='L')
@plac.flg("global_", "local alignment", abbrev='G')
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.opt("type_", "sequence type (nuc, pep)", choices=[DNA, PEP])
def run(match='', mismatch='', open_='', extend='', matrix='', local_=False, global_=False,
        semiglobal=False, type_='', vcf=False, table=False, vars=False, *sequences):
    # Select alignment mode
    if global_:
        mode = GLOBAL_ALIGN
    elif local_:
        mode = LOCAL_ALIGN
    elif semiglobal:
        mode = SEMIGLOBAL_ALIGN
    else:
        mode = SEMIGLOBAL_ALIGN

    # Select formatting mode
    if vcf:
        fmt = VCF_FMT
    elif table:
        fmt = TABLE_FMT
    elif vars:
        fmt = VARIANTS_FMT
    else:
        fmt = PAIRWISE_FMT

    # Input data sources.
    lines = []

    if not sys.stdin.isatty():
        lines.append(sys.stdin)

    # Command line will be second in line.
    lines.extend(sequences)

    # Names sequences that come from command line
    counter = cycle(string.ascii_uppercase)

    # Records to be aligned
    recs = []
    for text in lines:
        elems = parse_input(text, counter=counter)
        recs.extend(elems)

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
        if len(rec.seq) > MAXLEN:
            utils.error(f"Sequence {rec.id} is too long for this aligner: {len(rec)} > MAXLEN={MAXLEN:,}", stop=False)
            utils.error("We recommend that you use a different software.")

    # The first sequence is the target.
    target = recs[0]

    collect = []
    # Generate an alignment for each record
    for query in recs[1:]:

        par = models.Param(
            gap_open=open_,
            gap_extend=extend,
            matrix=matrix,
            match=match,
            mismatch=mismatch,
            type=type_,
            mode=mode,
            fmt=fmt,
            target=target,
            query=query
        )

        alns = align(target, query, par=par)

        for aln in alns:

            # This is how you unpack a BioPython pairwise result.
            seq1, trace, seq2 = format(aln).splitlines()

            # Wrap the resulting alignment into a sequence.
            target_aln = models.Sequence(title=target.name, seq=seq1)
            query_aln = models.Sequence(title=query.name, seq=seq2)

            # Pack all content into the representation.
            obj = models.Alignment(target=target_aln, query=query_aln, trace=trace, aln=aln, score=aln.score, par=par)

            # Collect all alignments into a datastructure
            collect.append(obj)

            # Keep the first alignment only
            break

    if par.fmt == VCF_FMT:
        models.format_vcf(collect)
    elif par.fmt == TABLE_FMT:
        models.format_table(collect)
    elif par.fmt == VARIANTS_FMT:
        models.format_variants(collect)
    else:
        models.format_pairwise(collect)


if __name__ == '__main__':
    # bio align AGATTACA GATCA
    run()
