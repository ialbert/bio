import io
import os
import string
import subprocess
import sys
import tempfile
from itertools import *
from subprocess import PIPE

import plac

from biorun import models

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import substitution_matrices
except ImportError as exc:
    print(f"# Error: {exc}", file=sys.stderr)
    print(f"# This program requires biopython", file=sys.stderr)
    print(f"# Install: conda install -y biopython>=1.79", file=sys.stderr)
    sys.exit(-1)

from biorun import utils

DNA, PEP = "DNA", "PEP"

TABLE_FMT, VCF_FMT, VARIANTS_FMT, PAIRWISE_FMT, FASTA_FMT = "table", "vcf", 'variants', 'pariwise', 'fasta'

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


def guess_type(text):
    if all_nuc(text):
        return NUCLEOTIDE
    elif all_pep(text):
        return PEPTIDE
    else:
        return None


def parse_input(handle, counter):
    if hasattr(handle, "read", ):
        # Object may already be a stream
        recs = SeqIO.parse(handle, "fasta")
    elif os.path.isfile(handle):
        # Perhaps it is a filename
        recs = SeqIO.parse(handle, "fasta")
    elif guess_type(handle):
        # Perhaps the sequence comes from command line
        name = next(counter)
        name = f"{name}"
        recs = [SeqRecord(id=name, name=name, description='', seq=Seq(handle))]
    else:
        utils.error(f"Invalid file/sequence: {handle}")
        recs = []

    return list(recs)


def safe_abs(value):
    try:
        return abs(float(value))
    except Exception as exc:
        utils.error(f"{exc} {value}")


def run_cmd(cmd):
    proc = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=PIPE)

    if proc.returncode != 0:
        print(proc.stdout.decode("UTF-8"))
        print(proc.stderr.decode("UTF-8"))
        print("-" * 10)
        print(cmd)
        sys.exit(1)
    return proc


def run_aligner(query, targets, par: models.Param):
    num1, fname1 = tempfile.mkstemp()
    num2, fname2 = tempfile.mkstemp()

    # Set alignment defaults
    par.gap_open = abs(par.gap_open)
    par.gap_extend = abs(par.gap_extend)

    mapper = {
        LOCAL_ALIGN: 'matcher',
        GLOBAL_ALIGN: 'stretcher',
        SEMIGLOBAL_ALIGN: 'needle',

    }

    program = mapper[par.mode]

    results = []
    try:
        path = matrix_path(par.matrix)
        SeqIO.write(query, fname1, format='fasta')
        for target in targets:
            SeqIO.write(target, fname2, format='fasta')

            cmd = f'{program} {fname1} {fname2} --filter -aformat3 fasta -gapopen {par.gap_open} -gapextend {par.gap_extend} -datafile {path}'
            proc = run_cmd(cmd)
            output = proc.stdout.decode("UTF-8").strip()
            stream = io.StringIO(output)
            recs = SeqIO.parse(stream, format='fasta')

            # Generate an alignment for the data
            query = next(recs)
            target = next(recs)

            aln = models.Alignment(query=query, target=target)

            # Read the errors
            text_err = proc.stderr.decode("UTF-8")

            if text_err:
                sys.stderr.write(text_err)

            results.append(aln)
    finally:
        os.remove(fname1)
        os.remove(fname2)

    return results


def matrix_path(matrix):
    """
    Attempts to figure out a matrix path
    """

    if matrix.upper() in ("DNA", "EDNAFULL"):
        matrix = "NUC.4.4"

    try:
        if os.path.isfile(matrix):
            path = matrix
        else:
            # Attempt to load the matrix as it generates a better error message
            mat = substitution_matrices.load(matrix)
            path = os.path.dirname(os.path.realpath(substitution_matrices.__file__))
            path = os.path.join(path, "data", matrix)

    except Exception as exc:
        valid = substitution_matrices.load()
        utils.error(f"error loading matrix: {matrix}", stop=False)
        #utils.error(f"{exc}", stop=False)
        utils.error(f"built-in matrices: {', '.join(valid)}")
        sys.exit()

    return path


def show_matrix(matrix):
    path = matrix_path(matrix)
    print(open(path).read(), end='')


@plac.pos("sequence", "sequences")
@plac.opt("open_", "gap open penalty", abbrev='o')
@plac.opt("extend", "gap extend penalty", abbrev='x')
@plac.opt("matrix", "matrix (default: NUC.4.4 for DNA, BLOSUM62 for PEP)", abbrev='m')
@plac.opt("output", "output vcf file", abbrev='O', choices=["pairwise", "vcf", "table", "vars", "fasta"])
@plac.flg("local_", "local alignment", abbrev='L')
@plac.flg("global_", "local alignment", abbrev='G')
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.opt("type_", "sequence type (nuc, pep)", choices=[DNA, PEP])
@plac.flg("fasta", "fasta output", abbrev='F')
@plac.flg("table", "tabular output", abbrev='T')
@plac.flg("vcf", "vcf output", abbrev='V')
@plac.flg("vars", "variant table output", abbrev='A')
@plac.flg("pairwise", "pairwise output", abbrev='P')
def run(open_=11, extend=1, matrix='', local_=False, global_=False,
        semiglobal=False, type_='', fasta=False, table=False, vcf=False, vars=False, pairwise=False, *sequences):

    for rec in sequences:
        if rec.startswith("-"):
            utils.error(f"unrecognized flag: {rec}")

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
    elif fasta:
        fmt = FASTA_FMT
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
    counter = cycle(string.ascii_lowercase)

    # Records to be aligned
    recs = []
    for text in lines:
        elems = parse_input(text, counter=counter)
        recs.extend(elems)

    # If only matrix is specified print it to the screen.
    if matrix and len(recs) == 0:
        show_matrix(matrix)
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

    # The first sequence is the queries
    query = recs[0]

    if type_:
        is_dna = type_ == DNA
    else:
        is_dna = all_nuc(query.seq)

    # Setting the default matrix
    if not matrix:
        matrix = 'NUC.4.4' if is_dna else "BLOSUM62"

    # Other sequences are targets.
    targets = recs[1:]

    par = models.Param(
        gap_open=open_,
        gap_extend=extend,
        matrix=matrix,
        type=type_,
        mode=mode,
        fmt=fmt,
        is_dna=is_dna
    )

    results = run_aligner(query=query, targets=targets, par=par)

    if par.fmt == VCF_FMT:
        models.format_vcf(results)
    elif par.fmt == TABLE_FMT:
        models.format_table(results)
    elif par.fmt == VARIANTS_FMT:
        models.format_variants(results)
    elif par.fmt == FASTA_FMT:
        models.format_fasta(results)
    else:
        models.format_pairwise(results, par=par)


if __name__ == '__main__':
    # bio align AGATTACA GATCA
    plac.call(run)
