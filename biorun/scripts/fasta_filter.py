"""
Implements utilities for filtering FASTA files.
"""
import sys, plac, operator

from Bio import SeqIO

def nointerrupt(func):
    """
    Intercept keyboard interrupts.
    """

    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except KeyboardInterrupt:
            sys.exit(0)

    return wrapper

def print_acc(records):
    for rec in records:
        print (rec.id)

@plac.opt('minL', "minimum lenght", type=int, abbrev="L")
@plac.opt('maxN', "maximum number of Ns", abbrev="N")
@plac.opt('maxX', "maximum number of Xs", abbrev="X")
@plac.flg('invert', "invert filtering action")
@plac.flg('acc', "print accession numbers only")
def main(minL=-1, maxN=-1, maxX=-1, invert=False, acc=False):
    # Parse the input string.
    stream = SeqIO.parse(sys.stdin, "fasta")

    # Apply the minimal lenght filter
    if minL > -1:
        op1 = operator.lt if invert else operator.ge
        stream = filter(lambda r: op1(len(r.seq), minL), stream)

    # Apply maxN filter.
    if maxN > -1:
        op2 = operator.gt if invert else operator.le
        stream = filter(lambda r: op2(r.seq.count("N"), maxN), stream)

    if maxX > -1:
        op2 = operator.gt if invert else operator.le
        stream = filter(lambda r: op2(r.seq.count("X"), maxX), stream)

    # Produce the output
    if acc:
        print_acc(stream)
    else:
        SeqIO.write(stream, sys.stdout, "fasta")


@nointerrupt
def run():
    """
    Entry point for the script.
    """
    plac.call(main)

if __name__ == '__main__':
    run()
