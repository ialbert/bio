"""
Generates tests from the test_builder.sh shell script.

Creates the test file called test_bio.py

Each line in the shell script will be a line in the

"""
from itertools import count
from textwrap import dedent
import os, sys, difflib
import plac
from biorun import main, const
from biorun.methods import align

# Test naming index.
counter = count(1)

# The path to the current file.
CURR_DIR = os.path.dirname(__file__)

# The default data directory.
DATA_DIR = os.path.join(CURR_DIR, "data")


def read(fname, datadir=DATA_DIR):
    """
    Reads a file in the datadir
    """
    path = os.path.join(datadir, fname) if datadir else fname
    text = open(path).read()
    return text


def run(cmd, capsys, out=None):
    """
    Runs a command and returns its out.
    """

    # Drop the leading command (bio)
    params = cmd.split()[1:]

    # Different functions to be called based on the command.
    if params and const.ALIGN in params:
        # Run the alignment tests.
        params.remove(const.ALIGN)
        assert plac.call(align.run, params) is None
    else:
        # Run converter commands.
        assert plac.call(main.converter, params) is None

    # Read the standard out
    stream = capsys.readouterr()
    result = stream.out

    # Check the output if we pass expected value here.
    if out:
        expect = read(out)
        if expect != result:
            lines1 = expect.splitlines()
            lines2 = result.splitlines()
            diffs = difflib.unified_diff(lines1, lines2)
            for diff in diffs:
                print(diff)
            assert result == expect

    return result


init = '''
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from generate import run
'''


def generate_tests(infile, outfile="test_bio.py"):
    """
    Generates tests from a shell script.
    """
    print(f"*** script {infile}")
    print(f"*** tests {outfile}")

    stream = open(infile)
    lines = map(lambda x: x.strip(), stream)
    lines = filter(lambda x: ">" in x, lines)
    lines = map(lambda x: tuple(x.split(">")), lines)
    lines = list(lines)

    collect = []
    for cmd, fname in lines:
        fname = fname.strip()
        patt = f"""
        def test_{next(counter)}(capsys):
            cmd = "{cmd}"
            run(cmd, capsys=capsys, out="{fname}")

        """
        collect.append(dedent(patt))

    fp = open(outfile, "wt")
    print(init, file=fp)
    print("".join(collect), file=fp)
    fp.close()

if __name__ == '__main__':
    infile = os.path.join(CURR_DIR, "test_bio_data.sh")
    outfile = os.path.join(CURR_DIR, "test_bio.py")
    generate_tests(infile=infile, outfile=outfile)
