import os, csv
from biorun import __main__ as bio
import plac
from biorun.align import pairwise

# The path to the current file.
__CURR_DIR = os.path.dirname(__file__)

__DATADIR = os.path.join(__CURR_DIR, "data")


def read_file(fname, datadir=__DATADIR):
    path = os.path.join(datadir, fname) if datadir else fname
    text = open(path).read()
    return text


def run_bio(cmd, capsys, output=None):
    """
    Runs a command and returns its output.
    """

    # Take the parameters only.
    params = cmd.split()[1:]

    # Run the command and assert its state.
    assert plac.call(bio.run, params) == None

    # Read the standard output
    stream = capsys.readouterr()
    result = stream.out

    if output:
        assert result[:1000] == output[:1000]

    return result


def test_fetch(capsys):
    cmd = "bio fetch NC_045512"
    run_bio(cmd, capsys=capsys)

def test_view(capsys):
    cmd = "bio view NC_045512"
    output = read_file("NC_045512.gb")
    run_bio(cmd, capsys=capsys, output=output)

def test_view_fasta(capsys):
    cmd = "bio view NC_045512 -f"
    output = read_file("NC_045512.fa")
    run_bio(cmd, capsys=capsys, output=output)


def main():
    pass


if __name__ == '__main__':
    main()
