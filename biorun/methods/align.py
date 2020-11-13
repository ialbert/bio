import warnings, sys, os
import biorun.libs.placlib as plac
from itertools import islice, count
import textwrap

from biorun import const
from biorun import utils, storage, objects
from biorun.models import jsonrec, fastarec

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

try:
    import parasail
except ImportError as exc:
    print(f"*** Warning: {exc}", file=sys.stderr)
    print(f"*** Please install parasail: conda install -y parasail-python", file=sys.stderr)
    sys.exit(-1)


# The default logging function.
logger = utils.logger


class Alignment():
    """
    A wrapper class to represent an alignment.
    """

    def __init__(self, query, target, trace,
                 gap_open=11, gap_extend=1, matrix='', mode='',
                 ichr='|', mchr='.', gchr=' ', schr=':', attrs={}):
        self.query = query
        self.target = target
        self.trace = trace
        self.len = len(trace)
        self.mode = mode
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.matrix = matrix
        self.len_query = self.len_ref = self.score = 0
        self.start_query = self.start_ref = self.end_query = self.end_ref = 0

        # Identity
        self.icount = self.trace.count(ichr)
        self.iperc = 100 * self.icount / self.len if self.len else 0

        # Similarity
        self.scount = self.icount + trace.count(schr)
        self.sperc = 100 * self.scount / self.len if self.len else 0

        # Mismatches
        self.mcount = trace.count(mchr)
        self.mperc = 100 * self.mcount / self.len if self.len else 0

        # Gaps
        self.gcount = trace.count(gchr)
        self.gperc = 100 * self.gcount / self.len if self.len else 0

        # Update with additional attributes.
        self.__dict__.update(attrs)

        self.counter = count(1)

    def print_wrapped(self, width=100, **kwargs):
        """
        Wraps and prints alignments
        """

        # Enforce a fixed width on each name.
        get = lambda name: f'{kwargs.get(name, ""):12.12s}'

        # Fetch the query names from additional attributes.
        q_name = get("q_name")
        p_name = get("p_name")
        t_name = get("t_name")

        header = f'''
        ### {next(self.counter)}: {q_name.strip()} vs {t_name.strip()} ###

        Length:\t{self.len} ({self.mode}) 
        Query:\t{self.len_query} [{self.start_query}, {self.end_query}]
        Target:\t{self.len_ref} [{self.start_ref}, {self.end_ref}]
        Score:\t{self.score}
        Ident:\t{self.icount}/{self.len} ({self.iperc:.1f}%)
        Simil:\t{self.scount}/{self.len} ({self.sperc:.1f}%)
        Gaps:\t{self.gcount}/{self.len} ({self.gperc:.1f}%)
        Matrix:\t{self.matrix}(-{self.gap_open}, -{self.gap_extend}) 
        '''

        header = textwrap.dedent(header)
        print(header)
        for start in range(0, len(self.trace), width):
            end = min(start + width, self.len)
            print(q_name, self.query[start:end])
            patt = f"{start + 1}"
            print(f"{patt:>12.12s} ", end="")
            print(self.trace[start:end], end="")
            print (f" {end}")
            print(t_name, self.target[start:end])
            print ("")


def get_matrix(seq, matrix):
    if not matrix:
        haspep = any(x for x in seq[:100] if x not in "ATGC")
        matrix = parasail.blosum62 if haspep else parasail.nuc44
    else:
        raise Exception("No matrix found")
    return matrix

# Maps constants to parasail functions.
FUNC_MAPPER = {
    # Smith-Waterman local alignment.
    const.LOCAL_ALIGN: parasail.sw_trace_scan,

    # Needleman-Wunsch global alignment.
    const.GLOBAL_ALIGN: parasail.nw_trace_scan,

    # Semiglobal alignment, no gap penalites on either end.
    const.SEMIGLOBAL_ALIGN: parasail.sg_trace_scan,
}

def parasail_align(qseq, tseq, param):
    q = str(qseq.seq)
    t = str(tseq.seq)

    # Guess matrix type
    matrix = get_matrix(t, param.matrix)

    mode = param.mode if (param.mode in FUNC_MAPPER) else const.SEMIGLOBAL_ALIGN

    # Select the function mapper
    func = FUNC_MAPPER[mode]

    res = func(q, t, param.gap_open, param.gap_extend, matrix=matrix)

    # Alignment must be traceback aware.
    t = res.traceback

    # Shortcuts to each field of the traceback.
    query, target, trace = t.query, t.ref, t.comp

    # Collect additional attributes
    attrs = dict()
    words = "score matches len_ref len_query cigar".split()
    for word in words:
        attrs[word] = getattr(res, word, '')

    # Populate the start coordinates
    coords = "start_ref end_ref start_query end_query".split()
    for coord in coords:
        attrs[coord] = getattr(res, coord, 0) + 1


    # Populate from cigar string.
    attrs['start_query'] = res.cigar.beg_query + 1
    attrs['start_ref'] = res.cigar.beg_ref + 1

    # Decode the CIGAR string
    attrs['cigar'] = res.cigar.decode.decode("ascii")

    # String name for the matrix
    mname = str(matrix.name.decode("ascii"))

    aln = Alignment(query=query, target=target, gap_open=param.gap_open, gap_extend=param.gap_extend,
                    trace=trace, attrs=attrs, matrix=mname, mode=param.mode)

    # For semiglobal alignment need to manually find the start/end from the pattern.
    aln.print_wrapped(q_name=qseq.id, t_name=tseq.id)


@plac.opt('start', "start coordinate ", type=int)
@plac.opt('end', "end coordinate")
@plac.opt('matrix', "scoring matrix", abbrev='M')
@plac.opt('gap_open', "scoring matrix", abbrev='o')
@plac.opt('gap_extend', "scoring matrix", abbrev='x')
@plac.flg('local_', "perform local alignment", abbrev='L')
@plac.flg('global_', "perform global alignment (zero end gap penalty)", abbrev='G')
@plac.flg('semiglobal', "perform a semiglobal alignment", abbrev='S')
@plac.flg('inter', "interactive mode, data from command line")
@plac.flg('protein', "use the translated protein sequences from the data")
@plac.flg('translate', "use the translated protein sequences from the data")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(start=1, end='', gap_open=11, gap_extend=1, local_=False, global_=False, semiglobal=False,
        protein=False, translate=False, inter=False, verbose=False, query='',  target=''):
    """
    Handles an alignment request.
    """

    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    # Ensure counter is reset.
    jsonrec.reset_counter()

    # Requires two inputs.
    if not (query and target):
        utils.error(f"Please specify both a QUERY and a TARGET")

    if global_:
        mode = const.GLOBAL_ALIGN
    elif local_:
        mode = const.LOCAL_ALIGN
    elif semiglobal:
        mode = const.SEMIGLOBAL_ALIGN
    else:
        mode = const.SEMIGLOBAL_ALIGN

    param1 = objects.Param(name=query, protein=protein, translate=translate,
                           start=start, end=end, gap_open=gap_open, gap_extend=gap_extend, mode=mode)
    param2 = objects.Param(name=target, protein=protein, translate=translate,
                           start=start, end=end, gap_open=gap_open, gap_extend=gap_extend, mode=mode)

    # Get the JSON data.
    param1.json = storage.get_json(param1.name, inter=inter, strict=True)
    param2.json = storage.get_json(param2.name, inter=inter, strict=True)

    for rec1 in param1.json:

        for rec2 in param2.json:

            qrecs = fastarec.get_fasta(rec1, param=param1)
            trecs = fastarec.get_fasta(rec2, param=param2)

            for qseq in qrecs:
                for tseq in trecs:
                    parasail_align(qseq=qseq, tseq=tseq, param=param1)


    # biopython_align(query=query, target=target, matrix=matrix)



def main():
    plac.call(run)


if __name__ == '__main__':
    main()
