import os, sys
import plac
import subprocess, difflib
from subprocess import PIPE
from tqdm import tqdm

join = os.path.join

# The path to the current file.
CURR_DIR = os.path.abspath(os.path.dirname(__file__))

# The default data directory.
DATA_DIR = os.path.join(CURR_DIR, "data")

# The run directory
RUN_NAME = "run"
RUN_DIR = os.path.join(CURR_DIR, RUN_NAME)

USAGE = join(CURR_DIR, "usage.sh")

def parse_commands(text, flag=False):
    lines = text.splitlines()
    lines = map(lambda x: x.strip(), lines)
    lines = filter(None, lines)
    lines = filter(lambda x: not x.startswith("#"), lines)
    if flag:
        lines = filter(lambda x: ">" in x, lines)
    lines = list(lines)
    return lines


def run(cmd):
    proc = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=PIPE)

    if proc.returncode != 0:
        print (proc.stdout.decode("UTF-8"))
        print(proc.stderr.decode("UTF-8"))
        print("-" * 10)
        print(cmd)

        sys.exit(1)
    return proc


INIT = f"""    
    ln -s ../data/alias.txt
    ln -s ../data/align_input.fa
    ln -s ../data/mafft.fa
"""


def run_setup():
    os.chdir(CURR_DIR)
    run(f"rm -rf {RUN_NAME}")
    run(f"mkdir {RUN_NAME}")
    os.chdir(RUN_DIR)
    init = parse_commands(INIT)
    for cmd in init:
        run(cmd)

def print_diff(expect, result):
    """
    User friendly diffs
    """
    lines1 = expect.splitlines()
    lines2 = result.splitlines()
    diffs = difflib.unified_diff(lines1, lines2)
    for diff in diffs:
        print(diff)

def main():

    run_setup()

    text = open(USAGE).read()
    cmds = parse_commands(text, flag=True)

    #cmds = cmds[:5]

    print(f"# Running {len(cmds)} tests:")

    cmds = tqdm(cmds)

    for cmd in cmds:
        fname = cmd.split(">")[1].strip()
        run(cmd)
        result = open(join(RUN_DIR, fname)).read()
        expect = open(join(DATA_DIR, fname)).read()
        if expect != result:
            print_diff(expect=expect, result=result)
            print (f"\n(cd data && {cmd})\n")
            sys.exit(1)

    print (f"# All tests completed")
if __name__ == '__main__':
    plac.call(main)
