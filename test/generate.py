"""
Generates tests from the test_builder.sh shell script.

Creates the test file called test_bio.py

Each line in the shell script will be a line in the

"""
from itertools import count
from textwrap import dedent
import os, sys, difflib
from biorun import main

# Test naming index.
COUNTER = count(1)

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


def run(cmd, capsys, fname=None):
    """
    Runs a command and returns its out.
    """

    # Override the system arguments.
    sys.argv = cmd.split()

    # Dispatch the commands.
    main.router()

    # Read the standard output.
    result = capsys.readouterr().out

    # Check the output if we pass expected value here.
    if fname:
        expect = read(fname)
        if expect != result:
            lines1 = expect.splitlines()
            lines2 = result.splitlines()
            diffs = difflib.unified_diff(lines1, lines2)
            print (cmd)
            print (f"File: {fname}")
            print ("-" * 10)
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
    # Test only bio commands.
    lines = filter(lambda x: x[:3] == "bio", lines)
    lines = list(lines)

    collect = []
    for line in lines:
        if ">" in line:
            cmd, fname = line.split(">")
            cmd = cmd.strip()
            fname = fname.strip()
            fname = f'"{fname}"'
        else:
            cmd, fname = line, None

        patt = f"""
        def test_{next(COUNTER)}(capsys):
            cmd = "{cmd}"
            run(cmd, capsys=capsys, fname={fname})
        """
        collect.append(dedent(patt))

    fp = open(outfile, "wt")
    print(init, file=fp)
    print("".join(collect), file=fp)
    fp.close()

if __name__ == '__main__':
    infile = os.path.join(CURR_DIR, "bio-examples.sh")
    outfile = os.path.join(CURR_DIR, "test_bio.py")
    generate_tests(infile=infile, outfile=outfile)
