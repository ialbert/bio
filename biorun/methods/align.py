import warnings, sys, os
import biorun.libs.placlib as plac
from itertools import islice, count
import textwrap

from biorun import const
from biorun import utils, storage, objects
from biorun.models import jsonrec, fastarec

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
except ImportError as exc:
    print(f"*** Please install biopython: conda install -y biopython", file=sys.stderr)
    print(f"*** Error: {exc}", file=sys.stderr)
    sys.exit(-1)

# The default logging function.
logger = utils.logger

class Alignment:
    """
    A wrapper class to represent an alignment.
    """

    def __init__(self,  qseq, tseq, aln, param):
        self.qseq = qseq
        self.tseq = tseq
        self.score = aln.score
        self.path = aln.path
        self.aln = aln
        self.param = param
        self.query = self.trace = self.target = ""
        self.len = self.len_query = self.len_target = 0
        self.start_query = self.start_target = self.end_query = self.end_target = 0
        self.icount = self.iperc = self.scount = self.sperc = 0
        self.mcount = self.mperc = self.gcount = self.gperc = 0

    def format_trace(self, ichr='|', mchr='.', gchr='-', schr=':'):

        # Format for the trace
        text = format(self.aln)

        # Extract the alignment elements.
        self.query, self.trace, self.target = text.splitlines()

        # Show only aligned regions for local and semiglobal alignments
        if self.param.mode in (const.LOCAL_ALIGN, const.SEMIGLOBAL_ALIGN):

            char = " " if self.param.mode == const.LOCAL_ALIGN else "-"
            lcount = len(self.trace) - len(self.trace.lstrip(char))
            rcount = len(self.trace) - len(self.trace.rstrip(char))

            if lcount:
                if self.query.startswith(char):
                    self.target = self.target[lcount:]
                else:
                    self.query = self.query[lcount:]

            if rcount:
                if self.query.endswith(char):
                    self.target = self.target[:-rcount]
                else:
                    self.query = self.query[:-rcount]

            self.target = self.target.strip(char)
            self.trace = self.trace.strip(char)
            self.query = self.query.strip(char)

        # Alignment length
        self.len = len(self.trace)

        # Identity
        self.icount = self.trace.count(ichr)
        self.iperc = 100 * self.icount / self.len if self.len else 0

        # Similarity
        self.scount = self.icount + self.trace.count(schr)
        self.sperc = 100 * self.scount / self.len if self.len else 0

        # Mismatches
        self.mcount = self.trace.count(mchr)
        self.mperc = 100 * self.mcount / self.len if self.len else 0

        # Gaps
        self.gcount = self.trace.count(gchr)
        self.gperc = 100 * self.gcount / self.len if self.len else 0

        # Unpack the paths
        t_path, q_path = self.aln.aligned

        self.t_start = t_path[0][0] + 1
        self.t_end = t_path[-1][-1] + 1

        self.q_start = q_path[0][0] + 1
        self.q_end = q_path[-1][-1] + 1


def biopython_align(qseq, tseq, param, table=False, strict=False):

    # Query and target sequences.
    q = str(qseq.seq)
    t = str(tseq.seq)

    aligner = Align.PairwiseAligner()

    if param.mode == const.LOCAL_ALIGN:
        aligner.mode = 'local'

    # Select the default substituion matrix.
    if not param.matrix:
        # Detect peptide sequences
        haspep = any(x for x in q[:100] if x not in "ATGC")
        param.matrix = 'BLOSUM62' if haspep else 'NUC.4.4'

    aligner.substitution_matrix = substitution_matrices.load(param.matrix)

    # Gap scoring.
    aligner.open_gap_score = -param.gap_open
    aligner.extend_gap_score = -param.gap_extend

    # End gap scoring.
    if strict:
        aligner.target_end_open_gap_score = -param.gap_open
        aligner.target_end_extend_gap_score = -param.gap_extend

        aligner.query_end_open_gap_score = -param.gap_open
        aligner.query_end_extend_gap_score = -param.gap_extend
    else:
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

    # Semiglobal will override strict mode.
    if param.mode == const.SEMIGLOBAL_ALIGN:
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

    alns = aligner.align(t, q)

    def builder(aln):
        rec = Alignment(qseq=qseq, tseq=tseq, aln=aln, param=param)
        return rec

    alns = map(builder, alns)

    for index, aln in enumerate(alns):
        if table:
            print_tabular(aln, index=index)
        else:
            print_pairwise(aln)


# Enforce a fixed width on each name.
def padded(value, right=False):
    value = str(value)
    return f'{value:>12.12s}' if right else f'{value:12.12s}'

def print_tabular(aln, index=0):
    # Format the trace.
    aln.format_trace()

    # Formatting the identity
    ident = f"{aln.icount}"
    iperc = f"{aln.iperc:.1f}"
    gperc = f"{aln.gperc:.1f}"
    mperc = f"{aln.mperc:.1f}"

    if index == 0:
        head = [
            "query", "target", "pident", "ident", "mism", "gaps", "score", "alen", "tlen", "tstart", "tend", "qlen", "qstart", "qend"
        ]
        print("\t".join(head))


    data = [
        aln.qseq.id, aln.tseq.id, iperc, aln.icount, aln.mcount, aln.gcount, aln.score,
            aln.len, len(aln.tseq), aln.t_start, aln.t_end, len(aln.qseq), aln.q_start, aln.q_end,
    ]
    data = map(str, data)
    print("\t".join(data))


def print_pairwise(aln, width=90):
    """
    A detailed visual pairwise alignment.
    """
    # Format the trace.
    aln.format_trace()

    # Formatting the identity
    ident = f"{aln.icount}({aln.iperc:.1f}%)"
    gaps = f"{aln.gcount}({aln.gperc:.1f}%)"
    mism = f"{aln.mcount}({aln.mperc:.1f}%)"
    alns  = f"Target={(aln.t_start,aln.t_end)}  Query={(aln.q_start,aln.q_end)}"

    header = f'''
    # Ident={ident}  Mis={mism}  Gaps={gaps}  {alns}  Length={aln.len:,}  Score={aln.score:0.1f}  {aln.param.matrix}({aln.param.gap_open},{aln.param.gap_extend})
    '''

    header = textwrap.dedent(header)

    query_id = padded(aln.qseq.id)
    target_id = padded(aln.tseq.id)

    print(header)
    for start in range(0, len(aln.trace), width):
        end = min(start + width, aln.len)
        trace_id = padded("")
        print(query_id, aln.query[start:end])
        print(f"{trace_id} {aln.trace[start:end]} {end}")
        print(target_id, aln.target[start:end])
        print("")


@plac.pos("query", "query sequence to align")
@plac.pos("target", "target sequence to align")
@plac.opt('start', "start coordinate ", type=int)
@plac.opt('end', "end coordinate")
@plac.opt('matrix', "scoring matrix", abbrev='M')
@plac.opt('gap_open', "gap open penalty", abbrev='o')
@plac.opt('gap_extend', "gap extend penalty", abbrev='x')
@plac.flg('local_', "perform local alignment", abbrev='L')
@plac.flg('global_', "perform global alignment (zero end gap penalty)", abbrev='G')
@plac.flg('semiglobal', "perform a semiglobal alignment", abbrev='S')
@plac.flg('inter', "interactive mode, data from command line")
@plac.flg('protein', "use the translated protein sequences from the data")
@plac.flg('translate', "translate the DNA into proteins")
@plac.flg('table', "generate an alignment table", abbrev='T')
@plac.flg('strict', "strict global alignment, apply end gap penalties", abbrev='R')
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(start=1, end='', gap_open=11, gap_extend=1, local_=False, global_=False, semiglobal=False,
        protein=False, translate=False, inter=False, table=False, strict=False, verbose=False, target=None, query=None):
    """
    Performs an alignment between the query and target.
    """

    # Alignments over this size will take a long time!
    MAX_LEN = 100000

    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    # Reset counter (needed for consistency during testing).
    jsonrec.reset_counter()

    # This method requires two inputs.
    if not (query and target):
        utils.error(f"Please specify both a QUERY and a TARGET")

    if global_:
        mode = const.GLOBAL_ALIGN
    elif local_:
        mode = const.LOCAL_ALIGN
    elif semiglobal:
        mode = const.SEMIGLOBAL_ALIGN
    else:
        mode = const.GLOBAL_ALIGN

    # A parameter for each record.
    param1 = objects.Param(acc=query, protein=protein, translate=translate,
                           start=start, end=end, gap_open=gap_open, gap_extend=gap_extend, mode=mode)

    param2 = objects.Param(acc=target, protein=protein, translate=translate,
                           start=start, end=end, gap_open=gap_open, gap_extend=gap_extend, mode=mode)

    # Get the JSON data.
    param1.json = storage.get_json(param1.acc, inter=inter, strict=True)
    param2.json = storage.get_json(param2.acc, inter=inter, strict=True)

    # Each data object may contain several records.
    for rec1 in param1.json:
        for rec2 in param2.json:

            qrecs = fastarec.get_fasta(rec1, param=param1)
            trecs = fastarec.get_fasta(rec2, param=param2)

            for qseq in qrecs:
                for tseq in trecs:

                    if (len(qseq) > MAX_LEN):
                        utils.error(f"query is longer than maximum: {len(qseq):,} > {MAX_LEN:,}")

                    if (len(tseq) > MAX_LEN):
                        utils.error(f"target sequence is longer than maximum: {len(tseq):,} > {MAX_LEN:,}")

                    biopython_align(qseq=qseq, tseq=tseq, param=param1, table=table, strict=strict)



def main():
    plac.call(run)


if __name__ == '__main__':
    main()
