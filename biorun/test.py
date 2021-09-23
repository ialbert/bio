import difflib
import os
import subprocess
import sys
from subprocess import PIPE

from tqdm import tqdm

from biorun.libs import placlib as plac


def join(*args):
    return os.path.abspath(os.path.join(*args))

# The path to the current file.
FILE_DIR = join(os.path.dirname(__file__))

# The default data directory.
DATA_DIR = join(FILE_DIR, "data")

RUN_DIR_NAME = "testrun"

# The run directory
RUN_DIR = join(os.path.expanduser("~"), ".bio")

os.makedirs(RUN_DIR, exist_ok=True)

USAGE = join(DATA_DIR, "usage.sh")

def parse_commands(text, flag=False):
    lines = text.splitlines()
    lines = map(lambda x: x.strip(), lines)
    lines = filter(None, lines)
    lines = filter(lambda x: not x.startswith("#"), lines)
    if flag:
        lines = filter(lambda x: ">" in x, lines)
    lines = list(lines)
    return lines


def shell(cmd):
    proc = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=PIPE)

    if proc.returncode != 0:
        print(proc.stdout.decode("UTF-8"))
        print(proc.stderr.decode("UTF-8"))
        print("-" * 10)
        print(cmd)

        sys.exit(1)
    return proc


INIT = f"""    
    ln -s '{FILE_DIR}/data/alias.txt'
    ln -s '{FILE_DIR}/data/align_input.fa'
    ln -s '{FILE_DIR}/data/mafft.fa'
    ln -s '{FILE_DIR}/data/file1.txt'
    ln -s '{FILE_DIR}/data/file2.txt'
"""


def test_setup():
    os.chdir(RUN_DIR)
    shell(f"rm -rf '{RUN_DIR_NAME}'")
    shell(f"mkdir '{RUN_DIR_NAME}'")
    os.chdir(RUN_DIR_NAME)

    init = parse_commands(INIT)
    for cmd in init:
        shell(cmd)


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
    test_setup()

    text = open(USAGE).read()
    cmds = parse_commands(text, flag=True)

    # cmds = cmds[:5]

    print(f"# Running {len(cmds)} tests:")

    cmds = tqdm(cmds)

    for cmd in cmds:
        fname = cmd.split(">")[1].strip()
        try:
            shell(cmd)
            result = open(join(RUN_DIR, RUN_DIR_NAME, fname)).read()
            expect = open(join(DATA_DIR, fname)).read()
        except Exception as exc:
            print(f"\n\n(cd {DATA_DIR} && {cmd})\n")
            print(f"*** error: {exc}")
            sys.exit(1)

        if expect != result:
            print(f"running: {cmd}")
            print_diff(expect=expect, result=result)
            print(f"\n\n(cd  {DATA_DIR} && {cmd})\n")
            sys.exit(1)

    print(f"# All tests completed")


def run():
    plac.call(main)


if __name__ == '__main__':
    run()
