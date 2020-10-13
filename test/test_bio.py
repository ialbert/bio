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
        if result != output:
            lines = result.splitlines()[:5]
            text = "\n".join(lines)
            print(text)
            assert False


    return result


def test_fetch(capsys):
    cmd = "bio fetch NC_045512"
    run_bio(cmd, capsys=capsys)

def test_view(capsys):
    cmd = "bio view NC_045512"
    output = read_file("NC_045512.gb")
    run_bio(cmd, capsys=capsys, output=output)

def test_view_list(capsys):
    cmd = "bio list"
    run_bio(cmd, capsys=capsys)

def test_view_fasta(capsys):
    cmd = "bio view NC_045512 --fasta"
    output = read_file("NC_045512.fa")
    run_bio(cmd, capsys=capsys, output=output)

def test_view_fasta_start(capsys):
    cmd = "bio view NC_045512 --fasta --rename foo --start 10 --end 20"
    output = read_file("parts/fasta-start.fa")
    run_bio(cmd, capsys=capsys, output=output)

def test_view_fasta_type(capsys):
    cmd = "bio view NC_045512 --fasta --type CDS"
    output = read_file("parts/CDS.fa")
    run_bio(cmd, capsys=capsys, output=output)

def test_view_fasta_type_start(capsys):
    cmd = "bio view NC_045512 --fasta --type gene --end 10"
    output = read_file("parts/gene-start.fa")
    run_bio(cmd, capsys=capsys, output=output)

def test_view_gff1(capsys):
    cmd = "bio view NC_045512 --gff"
    output = read_file("NC_045512.gff")
    run_bio(cmd, capsys=capsys, output=output)


def test_view_gff_name(capsys):
    cmd = "bio view NC_045512 --gff --gene S"
    output = read_file("parts/gene.gff")
    run_bio(cmd, capsys=capsys, output=output)

def test_view_gff_start(capsys):
    cmd = "bio view NC_045512 --gff  --start 21563 --end 21565"
    output = read_file("parts/start.gff")
    run_bio(cmd, capsys=capsys, output=output)


def test_view_gff_type(capsys):
    cmd = "bio view NC_045512 --gff  --type CDS"
    output = read_file("parts/type.gff")
    run_bio(cmd, capsys=capsys, output=output)


def main():
    pass


if __name__ == '__main__':
    main()
