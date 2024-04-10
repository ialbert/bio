import difflib
import os
import shutil
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

RUN_DIR_NAME = "test"

# The run directory
RUN_DIR = join(os.path.expanduser("~"), ".bio", RUN_DIR_NAME)

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


def test_setup():

    names = [
        'genomes.gb', 'alias.txt', 'align_input.fa', 'mafft.fa', 'file1.txt', 'file2.txt'
    ]
    for name in names:
        src = f'{FILE_DIR}/data/{name}'
        dest = f'{RUN_DIR}/{name}'
        shutil.copyfile(src, dest)

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

    os.chdir(RUN_DIR)

    print(f"# Running {len(cmds)} tests:")

    cmds = tqdm(cmds)

    # Run the commands.
    for cmd in cmds:
        fname = cmd.split(">")[1].strip()
        try:
            shell(cmd)
            result = open(join(RUN_DIR, fname)).read()
            expect = open(join(DATA_DIR, fname)).read()
        except Exception as exc:
            print(f"\n\n(cd {DATA_DIR} && {cmd})\n")
            print(f"*** error: {exc}")
            sys.exit(1)

        # Initialize the SequenceMatcher with the file contents
        matcher = difflib.SequenceMatcher(None, expect, result)

        # Publihsed data may be slightly different.
        ratio = matcher.ratio()

        if ratio < 0.9:
            print(f"running: {cmd}")
            print_diff(expect=expect, result=result)
            print(f"\n\n(cd  {DATA_DIR} && {cmd})\n")
            sys.exit(1)

    print(f"# All tests completed")


def run():
    plac.call(main)


if __name__ == '__main__':
    run()
