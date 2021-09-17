#
# Sanity check to verify that the alignments work
#
import difflib
import os
import random
import subprocess
import sys
from subprocess import PIPE

join = os.path.join

# The path to the current file.
CURR_DIR = os.path.abspath(os.path.dirname(__file__))

# The default data directory.
DATA_DIR = os.path.join(CURR_DIR, "data")

# The run directory
RUN_NAME = "run"
RUN_DIR = os.path.join(CURR_DIR, RUN_NAME)


def run(cmd):
    proc = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=PIPE)

    if proc.returncode != 0:
        print(proc.stdout.decode("UTF-8"))
        print(proc.stderr.decode("UTF-8"))
        print("-" * 10)
        print(cmd)

        sys.exit(1)
    return proc


def print_diff(expect, result):
    """
    User friendly diffs
    """
    lines1 = expect.splitlines()
    lines2 = result.splitlines()
    diffs = difflib.unified_diff(lines1, lines2)
    for diff in diffs:
        print(diff)

LEN = 100
MUT = 20

DNA = "ATGC"
PEP = "ARNDBCEQZGHILKMFPSTWYX"
SEQ = DNA

def main():
    from random import choice

    target = [choice(SEQ) for x in range(LEN)]
    query = list(target)
    idx = [ random.randint(0, LEN-1) for x in range(MUT)]
    for i in idx:
        query[i] = random.choice(SEQ)

    target = "".join(target)
    query = "".join(query)

    inp1 = join(RUN_DIR, "in1.fa")
    inp2 = join(RUN_DIR, "in2.fa")

    fp1 = open(inp1, 'wt')
    fp1.write(f">a\n{query}")
    fp1.close()

    fp2 = open(inp2, 'wt')
    fp2.write(f">b\n{target}")
    fp2.close()

    out1 = join(RUN_DIR, "aln1.fa")
    out2 = join(RUN_DIR, "aln2.fa")

    cmd1 = f'needle {inp1} {inp2} -gapopen 11 -gapextend 1 -filter -aformat3 fasta > {out1}'

    cmd2 = f"bio align {inp1} {inp2} -open 11 -extend 1 --fasta > {out2}"

    proc = run(cmd1)
    proc = run(cmd2)

    result = open(join(RUN_DIR, out1)).read()

    expect = open(join(DATA_DIR, out2)).read()

    # Skip cases with multiople valid alignments.
    if "WARNING" in expect:
        return

    if expect != result:
        #print_diff(expect=expect, result=result)

        c1 = cmd1.split('>')[0]
        c2 = cmd2.split('>')[0]

        print ("")
        print (query)
        print(target)
        print ("")

        print(f"{c1}")

        print(f"{c2}")

        print ("")
        sys.exit(1)


if __name__ == '__main__':
    from tqdm import tqdm
    steps = range(100)
    steps = tqdm(steps)
    for i in steps:
        main()
