
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from generate import run


def test_1(capsys):
    cmd = "bio ncov --delete"
    run(cmd, capsys=capsys, out=None)

def test_2(capsys):
    cmd = "bio NC_045512 --fetch --rename ncov --seqid ncov"
    run(cmd, capsys=capsys, out=None)

def test_3(capsys):
    cmd = "bio ncov"
    run(cmd, capsys=capsys, out="ncov.json")

def test_4(capsys):
    cmd = "bio ncov --fasta"
    run(cmd, capsys=capsys, out="ncov.fa")

def test_5(capsys):
    cmd = "bio ncov --gff"
    run(cmd, capsys=capsys, out="ncov.gff")

def test_6(capsys):
    cmd = "bio ncov --match ORF1ab --type gene"
    run(cmd, capsys=capsys, out="match.json")

def test_7(capsys):
    cmd = "bio ncov --gff --gene S"
    run(cmd, capsys=capsys, out="gene1.gff")

def test_8(capsys):
    cmd = "bio ncov --gff  --start 10000 --end 20000"
    run(cmd, capsys=capsys, out="overlap.gff")

def test_9(capsys):
    cmd = "bio ncov --gff  --type CDS"
    run(cmd, capsys=capsys, out="type.gff")

def test_10(capsys):
    cmd = "bio ncov --fasta --seqid foo --start 10 --end 20"
    run(cmd, capsys=capsys, out="fasta-start.fa")

def test_11(capsys):
    cmd = "bio ncov --fasta --type CDS"
    run(cmd, capsys=capsys, out="CDS.fa")

def test_12(capsys):
    cmd = "bio ncov --fasta --type gene --end 10"
    run(cmd, capsys=capsys, out="gene-start.fa")

def test_13(capsys):
    cmd = "bio ncov --protein --start -10"
    run(cmd, capsys=capsys, out="protein-end.fa")

def test_14(capsys):
    cmd = "bio ncov --translate --type CDS"
    run(cmd, capsys=capsys, out="translate.fa")

def test_15(capsys):
    cmd = "bio ncov:S --fasta --protein --seqid foo"
    run(cmd, capsys=capsys, out="s_prot_foo.fa")

def test_16(capsys):
    cmd = "bio MN996532 --fetch --rename ratg13 --seqid ratg13"
    run(cmd, capsys=capsys, out=None)

def test_17(capsys):
    cmd = "bio align ncov ratg13 --end 200"
    run(cmd, capsys=capsys, out="align-dna.txt")

def test_18(capsys):
    cmd = "bio align ncov:S ratg13:S --end 80"
    run(cmd, capsys=capsys, out="align-dna-s.txt")

def test_19(capsys):
    cmd = "bio align ncov:S ratg13:S --protein"
    run(cmd, capsys=capsys, out="align-protein-s.txt")

def test_20(capsys):
    cmd = "bio align ncov:S ratg13:S --end 80 --translate"
    run(cmd, capsys=capsys, out="align-translated-s.txt")

def test_21(capsys):
    cmd = "bio align THISLINE ISALIGNED  -i --local"
    run(cmd, capsys=capsys, out="align-local.txt")

def test_22(capsys):
    cmd = "bio align THISLINE ISALIGNED -i --global"
    run(cmd, capsys=capsys, out="align-global.txt")

def test_23(capsys):
    cmd = "bio align THISLINE ISALIGNED -i --semiglobal"
    run(cmd, capsys=capsys, out="align-semiglobal.txt")

