import os
import sys
from itertools import *

from biorun import parser, convert
from biorun.libs import placlib as plac
from . import models

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import PairwiseAligner
    from Bio.Align import substitution_matrices
except ImportError as exc:
    print(f"# Error: {exc}", file=sys.stderr)
    print(f"# This program requires biopython", file=sys.stderr)
    print(f"# Install: conda install -y biopython>=1.79", file=sys.stderr)
    sys.exit(-1)

from biorun import utils

TABLE_FMT, VCF_FMT, VARIANTS_FMT, PAIRWISE_FMT = "table", "vcf", 'variants', 'pariwise'

LOCAL_ALIGN, GLOBAL_ALIGN, SEMIGLOBAL_ALIGN = "local", "global", "semiglobal"

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


def safe_abs(value):
    try:
        return abs(float(value))
    except Exception as exc:
        utils.error(f"{exc} {value}")


def align(target, query, par: models.Param):
    # Query and target sequences.
    t = str(target.seq).upper()
    q = str(query.seq).upper()

    aligner = PairwiseAligner()

    # Select local mode. Global, semiglobal are about scoring.
    if par.mode == LOCAL_ALIGN:
        aligner.mode = 'local'
    else:
        aligner.mode = 'global'

    # Recognize a few shorter matrix names.
    par.matrix = 'NUC.4.4' if par.matrix in ('DNA', 'EDNAFULL') else par.matrix
    par.matrix = 'BLOSUM62' if par.matrix in ('PROT', 'PEP') else par.matrix

    # Automatic detection of nucleotide vs peptide
    par.is_dna = all_nuc(t)

    # Peptides get the BLOSUM matrix as default.
    if not par.matrix and not par.is_dna:
        par.matrix = "BLOSUM62"

    if par.matrix:
        # Load the substitutions matrix
        aligner.substitution_matrix = get_matrix(par.matrix)
    else:
        # Nucleotide defaults.
        aligner.match_score = safe_abs(par.match)
        aligner.mismatch_score = -safe_abs(par.mismatch)

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


@plac.pos("sequence", "sequences")
@plac.opt("open_", "gap open penalty", type=int, abbrev='o')
@plac.opt("extend", "gap extend penalty", type=float, abbrev='x')
@plac.opt("match", "alignment match (DNA only)", abbrev='m')
@plac.opt("mismatch", "alignment mismatch (DNA only)", abbrev='i')
@plac.opt("matrix", "matrix default: NUC4.4. or BLOSUM62)", abbrev='M')
@plac.flg("vcf", "output vcf file", abbrev='V')
@plac.flg("table", "output in tabular format", abbrev="T")
@plac.flg("diff", "output mutations", abbrev="d")
@plac.flg("pile", "output pileup", abbrev="p")
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.flg("fasta", "output variant columns", abbrev="F")
@plac.flg("local_", "local alignment", abbrev='L')
@plac.flg("global_", "local alignment", abbrev='G')
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.flg("all_", "show all alignments", abbrev='A')
def run(open_=11, extend=1, matrix='', local_=False, global_=False, match=1, mismatch=2,
        semiglobal=False, vcf=False, table=False, diff=False, pile=False, fasta=False, all_=False, *sequences):
    # Select alignment mode
    if global_:
        mode = GLOBAL_ALIGN
    elif local_:
        mode = LOCAL_ALIGN
    elif semiglobal:
        mode = SEMIGLOBAL_ALIGN
    else:
        mode = SEMIGLOBAL_ALIGN

    recs = parser.get_records(sequences)

    # Keep only source annotated sequences (genome level)
    recs = filter(convert.source_only, recs)

    # Get the sequence records.
    recs = list(recs)

    # If only matrix is specified print it to the screen.
    if matrix and len(recs) == 0:
        get_matrix(matrix, show=True)
        sys.exit()

    # Sequences must be present to be aligned.
    if len(recs) < 2:
        utils.error("Need at least two sequences to align")

    # The first sequence is the query
    target = recs[0]

    # Resulting alignments
    collect = []

    par = models.Param(
        gap_open=open_,
        gap_extend=extend,
        matrix=matrix,
        mode=mode,
        match=match,
        mismatch=mismatch,
    )

    # Keeping people from accidentally running alignments that are too large.
    MAXLEN = 50000

    # Generate an alignment for each target
    for query in recs[1:]:

        if (len(query.seq) > MAXLEN) and (len(target.seq) > MAXLEN):
            utils.error(f"Both sequences {query.id} and {target.id} appear to be longer than {MAXLEN:,} bases.", stop=False)
            utils.error("We recommend that you use a different software.")

        alns = align(query, target, par=par)

        alns = islice(alns, 1000)

        for aln in alns:

            # Unpack a BioPython pairwise result.
            seq1, seq2 = Seq(aln[0]), Seq(aln[1])

            # Wrap the resulting alignment into a sequence.
            query_aln = SeqRecord(id=query.id, name=query.name, description='', seq=seq1)
            target_aln = SeqRecord(id=target.id, name=target.name, description='', seq=seq2)

            # Pack all content into the representation.
            obj = models.Alignment(query=query_aln, target=target_aln, score=aln.score)

            # Collect all alignments into a datastructure
            collect.append(obj)

            # Keep the first alignment only
            if not all_:
                break

        count = sum(1 for x in alns)

    # Select formatting mode
    if vcf:
        models.format_vcf(collect)
    elif table:
        models.format_table(collect)
    elif diff:
        models.format_mutations(collect)
    elif fasta:
        models.format_fasta(collect)
    elif pile:
        models.format_pile(collect)
    else:
        models.format_pairwise(collect, par=par)

    if count > 1 and not all_:
        print(f"# {count + 1} identically scoring alignments! Pass the flag -A to see them all", file=sys.stderr)


if __name__ == '__main__':
    # Identically scoring alignments
    # bio align TGCGGGGGAAA GACGTGTGTCGTG - A - -fasta

    # bio align AGATTACA GATCA
    run()
