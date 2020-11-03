import os, csv
from biorun import main
import plac
from biorun.align import pairwise
import pytest

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
        assert result == out
        # Print only a subsection of the file
        #if result != out:
        #    lines = result.splitlines()[:50]
        #    text = "\n".join(lines)
        #    print(text)
        #    assert False

    return result


def test_delete(capsys):
    cmd = "bio SARS2 --delete"
    run(cmd, capsys=capsys)


def test_empty(capsys):
    cmd = "bio SARS2"
    with pytest.raises(SystemExit) as e:
        run(cmd, capsys=capsys)
    assert e.type == SystemExit
    assert e.value.code == 1

def test_fetch(capsys):
    cmd = "bio NC_045512 --fetch --rename SARS2 --seqid SARS2"
    run(cmd, capsys=capsys)

def test_list(capsys):
    cmd = "bio --list"
    run(cmd, capsys=capsys)

def test_view(capsys):
    cmd = "bio SARS2"
    out = read("SARS2.json")
    run(cmd, capsys=capsys, out=out)

def test_view_match(capsys):
    cmd = "bio SARS2 --match ORF1ab --type gene "
    out = read("parts/match.json")
    run(cmd, capsys=capsys, out=out)

def test_view_fasta(capsys):
    cmd = "bio SARS2 --fasta"
    out = read("SARS2.fa")
    run(cmd, capsys=capsys, out=out)

def test_view_fasta_start(capsys):
    cmd = "bio SARS2 --fasta --seqid foo --start 10 --end 20"
    out = read("parts/fasta-start.fa")
    run(cmd, capsys=capsys, out=out)

def test_protein_end(capsys):
    cmd = "bio SARS2 --protein --start -10"
    out = read("parts/protein-end.fa")
    run(cmd, capsys=capsys, out=out)

def test_view_fasta_type(capsys):
    cmd = "bio SARS2 --fasta --type CDS"
    out = read("parts/CDS.fa")
    run(cmd, capsys=capsys, out=out)

def test_view_fasta_type_start(capsys):
    cmd = "bio SARS2 --fasta --type gene --end 10"
    out = read("parts/gene-start.fa")
    run(cmd, capsys=capsys, out=out)

def test_view_gff1(capsys):
    cmd = "bio SARS2 --gff"
    out = read("SARS2.gff")
    run(cmd, capsys=capsys, out=out)

def test_view_gff_name(capsys):
    cmd = "bio SARS2 --gff --gene S"
    out = read("parts/gene.gff")
    run(cmd, capsys=capsys, out=out)

def test_view_gff_start(capsys):
    cmd = "bio SARS2 --gff  --start 10000 --end 20000"
    out = read("parts/overlap.gff")
    run(cmd, capsys=capsys, out=out)

def test_view_gff_type(capsys):
    cmd = "bio SARS2 --gff  --type CDS"
    out = read("parts/type.gff")
    run(cmd, capsys=capsys, out=out)

if __name__ == '__main__':
    pass
