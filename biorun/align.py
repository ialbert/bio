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

TABLE_FMT, VCF_FMT, VARIANTS_FMT = "table", "vcf", 'variants'

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

    print(f"# {label}: {fmt.target.name} ({len(fmt.target.seq):,}) vs {fmt.query.name} ({len(fmt.query.seq):,}) score = {fmt.score}")
    print(f"# Alignment: pident={fmt.pident:0.1f}% len={fmt.tlen} ident={fmt.ident} mis={fmt.mis} del={fmt.dels} ins={fmt.ins}")

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


def table_fmt(fmt, sep="\t"):
    data = [
        f"{fmt.target.name}", f"{fmt.query.name}",
        f"{fmt.score}", f"{fmt.pident:0.1f}", f"{fmt.tlen}",
        f"{fmt.ident}", f"{fmt.mis}", f"{fmt.dels}", f"{fmt.ins}"
    ]
    line = sep.join(data)
    print(line)


MATCH, SNP, INS, DEL = 'M', 'SNP', 'INS', 'DEL'


def find_variants(ref, tgt):
    """
    Takes two aligned sequences and tabulates what variants can be found at a given position.
    This is the trickiest to get right ... mostly works but may still have bugs.

    Returns a dictionary keyed by position where each position describes the  variant as a tuple.
    """

    # Peptides are named differently
    is_pep = guess_type(ref.seq) == PEPTIDE


    stream = zip(count(), ref.seq, tgt.seq)

    # stream = islice(stream, 50000)

    collect = []
    variants = []
    lastop = None
    pos = 0
    for idx, a, b in stream:

        # Multiple sequence alignment
        if a == '-' and b == '-':
            continue

        if a == b:
            # Matches
            pos += 1
            op = MATCH

        elif b == '-':
            # Deletion from query.
            pos += 1
            op = DEL

        elif a == '-':
            # Insertion into query
            op = INS

        elif a != b:
            # Mismatching bases.
            pos += 1
            op = SNP

        else:
            raise Exception(f"Should never hit this (sanity check): {a} vs {b}")

        # Collect variants when the operator changes.
        if lastop != op:
            if collect:
                if lastop != MATCH:
                    variants.append((lastop, collect))
            # Reset the collector
            collect = []

        # Collect the positions that have been visited
        collect.append((idx, pos, a, b))

        # Keep track of the last previous operation
        lastop = op

    # Collect last element
    if collect:
        variants.append((lastop, collect))

    # This is necessary because consecutive variants may overlap SNP + INSERT for example
    vcfdict = dict()

    for key, elems in variants:

        if key == MATCH:
            continue

        # Lenght of variant
        size = len(elems)

        idx, pos = elems[0][0], elems[0][1]

        # Query starts with insertion
        pos = 1 if pos < 1 else pos

        base = ref.seq[idx:idx + size].strip('-')
        alt = tgt.seq[idx:idx + size].strip('-')

        info = f"TYPE={key}"

        if key == SNP:
            # Mismatches printed consecutively
            for idx, pos, base, alt in elems:
                name = f"{base}{pos}{alt}"
                value = [ref.name, str(pos), name, base, alt, ".", "PASS", info, "GT", "1"]
                vcfdict[pos] = value


        elif key == INS or key == DEL:
            # Handle insertions and deletions.

            # Push back on POS if it is not 1.
            if idx > 0:
                # For insertions the pos does not advance
                if key == DEL:
                    pos = pos - 1
                base = ref.seq[idx - 1] + base
                alt = tgt.seq[idx - 1] + alt
            else:
                # When there is no preceding base the last base must be used
                lastidx = elems[-1][0] + 1
                base = base + ref.seq[lastidx]
                alt = alt + tgt.seq[lastidx]

            alt = alt or '.'
            base = base or '.'

            if key == DEL:
                name = f"{pos}del{base[1:]}"
            else:
                name = f"{pos}ins{alt[1:]}"

            value = [ref.name, str(pos), name, base, alt, ".", "PASS", info, "GT", "1"]
            vcfdict[pos] = value

    return vcfdict

def variants_fmt(fmt):
    vcfdict = find_variants(fmt.seqA, fmt.seqB)
    for value in vcfdict.values():
        name = value[2]
        pos = value[1]
        base = value[3]
        alt = value[4]
        size = '0'
        info = value[7]
        if "SNP" in info:
            info = 'snp'
            size = 1
        elif 'DEL' in info:
            info = 'del'
            base = base[1:] if pos != '1' else base[:1]
            alt = '-' * len(base)
            size = len(base)

        elif 'INS' in info:
            info = 'ins'
            alt = alt[1:] if pos != '1' else alt[:1]
            base = '-' * len(alt)
            size = len(alt)

        data = [ pos, info,  base, alt, str(size) ]
        print("\t".join(data))

def vcf_fmt(fmt):

    vcfdict = find_variants(fmt.seqA, fmt.seqB)

    ref, query = fmt.seqA, fmt.seqB

    print('##fileformat=VCFv4.2')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FILTER=<ID=PASS,Description="All filters passed">')
    print('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">')
    print(f'##contig=<ID={ref.name},length={len(ref.seq.strip("-"))},assembly={ref.name}>')
    print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{query.name}")

    for value in vcfdict.values():
        print("\t".join(value))


def align(target, query, par):
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
        print("\t".join(elems))


@plac.pos("sequence", "sequences")
@plac.opt("match", "match", type=int, abbrev='m')
@plac.opt("mismatch", "mismatch", type=int, abbrev='s')
@plac.opt("open_", "gap open penalty", type=int, abbrev='o')
@plac.opt("extend", "gap extend penalty", type=int, abbrev='x')
@plac.opt("matrix", "matrix", abbrev='M')
@plac.flg("vcf", "output vcf file", abbrev='V')
@plac.flg("table", "output as a table", abbrev="T")
@plac.flg("variants", "output as variants", abbrev="A")
@plac.flg("local_", "local alignment", abbrev='L')
@plac.flg("global_", "local alignment", abbrev='G')
@plac.flg("semiglobal", "local alignment", abbrev='S')
@plac.opt("type_", "sequence type (nuc, pep)", choices=[DNA, PEP])
def run(open_=6, extend=1, matrix='', match=1, mismatch=2, local_=False, global_=False,
        semiglobal=False, type_='', vcf=False, table=False, variants=False, *sequences):
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

    if vcf:
        par.format = VCF_FMT
    elif table:
        par.format = TABLE_FMT
    elif variants:
        par.format = VARIANTS_FMT
    else:
        par.format = pairwise_fmt

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

    # This names sequences that come from command line
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
            utils.error("We recommend that you use a different software.", stop=False)
            utils.error(f"Sequence {rec.id} is too long for this aligner: {len(rec)} > MAXLEN={MAXLEN:,}")

    target = recs[0]

    # Select result formatter
    if par.format == TABLE_FMT:
        header = "target query score len pident match mism ins del"
        formatter = table_fmt
    elif par.format == VCF_FMT:
        header = ''
        formatter = vcf_fmt
    elif par.format == VARIANTS_FMT:
        header = "pos type target query len"
        formatter = variants_fmt

    else:
        header = ''
        formatter = pairwise_fmt

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
