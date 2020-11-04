"""
Generates tests from the builder script shell script.

Creates the test file called test_bio.py

"""
from itertools import count
from textwrap import dedent
counter = count(1)

init = '''
#
# This file was generated automatically!
#

import os
import plac
from biorun import main 

# The path to the current file.
__CURR_DIR = os.path.dirname(__file__)

__DATADIR = os.path.join(__CURR_DIR, "data")

def read(fname, datadir=__DATADIR):
    path = os.path.join(datadir, fname) if datadir else fname
    text = open(path).read()
    return text

def run(cmd, capsys, out=None):
    """
    Runs a command and returns its out.
    """

    # Take the parameters only.
    params = cmd.split()[1:]

    # Run the command and assert its state.
    assert plac.call(main.base_runner, params) == None

    # Read the standard out
    stream = capsys.readouterr()
    result = stream.out

    if out:
        assert result == read(out)

    return result
'''

def run():
    stream = open("data/build_data.sh")
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

    fp = open("test_bio.py", "wt")
    print (init, file=fp)
    print ("".join(collect), file=fp)
    fp.close()

if __name__ == '__main__':
    run()