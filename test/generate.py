"""
Generates tests from the test_builder.sh shell script.

Creates the test file called test_bio.py

Each line in the shell script will be a line in the

"""
from itertools import count
from textwrap import dedent
import os, sys, difflib

# Test naming index.
COUNTER = count(1)

# The path to the current file.
CURR_DIR = os.path.dirname(__file__)

# The default data directory.
DATA_DIR = os.path.join(CURR_DIR, "data")

# The run directory
RUN_NAME = "run"
RUN_DIR = os.path.join(CURR_DIR, RUN_NAME)

def init_dirs():
    """
    Initializes the test run
    """
    # Change to current directory.        
    os.chdir(CURR_DIR)

    # Clean up the test directory.
    cmd = f"rm -rf {RUN_NAME}"
    os.system(cmd)

    # Make the run directory
    cmd = f"mkdir -p {RUN_NAME}"
    os.system(cmd)

    # Switch to run directory
    os.chdir(RUN_DIR)

join = os.path.join

def diff(cmd, expect, result):
    """
    User friendly diffs
    """
    lines1 = expect.splitlines()
    lines2 = result.splitlines()
    diffs = difflib.unified_diff(lines1, lines2)
    print (cmd)
    print ("-" * 10)
    for diff in diffs:
        print(diff)

def run(cmd):
    """
    Runs a command and checks the output if it has one.
    """
    # Run the command
    os.system(cmd)

    # If it produces output verify the correctness
    parts = cmd.split(">")
    if len(parts) == 2:
        fname = parts[-1].strip()
        result = open(join(RUN_DIR, fname)).read() 
        expect = open(join(DATA_DIR, fname)).read()
        if expect != result:
            diff(cmd=cmd, expect=expect, result=result)
            assert result == expect

init = '''
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from generate import *

# Initialize directories.
init_dirs()
'''


def generate_tests(infile, outfile="test_bio.py"):
    """
    Generates tests from a shell script.
    """
    print(f"*** script {infile}")
    print(f"*** tests {outfile}")

    collect = []
    stream = open(infile)
    stream = filter(lambda x: not x.startswith('#'), stream)
    stream = filter(lambda x: x.strip(), stream)
    for cmd in stream:
        cmd = cmd.strip()
        patt = f"""
        def test_{next(COUNTER)}(capsys):
            run("{cmd}")
        """
        collect.append(dedent(patt))

    fp = open(outfile, "wt")
    print(init, file=fp)
    print("".join(collect), file=fp)
    fp.close()

if __name__ == '__main__':

    infile = os.path.join(CURR_DIR, "all_tests.sh")    
    outfile = os.path.join(CURR_DIR, "test_bio.py")

    # Generate tests from a subset of
    if len(sys.argv) > 1:
        infile = sys.argv[1]

    generate_tests(infile=infile, outfile=outfile)
