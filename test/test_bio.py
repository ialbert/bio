
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
    cmd = "bio ncov --genbank"
    run(cmd, capsys=capsys, out="ncov.gb")

def test_5(capsys):
    cmd = "bio ncov --fasta"
    run(cmd, capsys=capsys, out="ncov.fa")

def test_6(capsys):
    cmd = "bio ncov --gff"
    run(cmd, capsys=capsys, out="ncov.gff")

def test_7(capsys):
    cmd = "bio ncov --match ORF1ab --type gene"
    run(cmd, capsys=capsys, out="match.json")

def test_8(capsys):
    cmd = "bio ncov --gff --gene S"
    run(cmd, capsys=capsys, out="gene1.gff")

def test_9(capsys):
    cmd = "bio ncov --gff  --start 10000 --end 20000"
    run(cmd, capsys=capsys, out="overlap.gff")

def test_10(capsys):
    cmd = "bio ncov --gff  --type CDS"
    run(cmd, capsys=capsys, out="type.gff")

def test_11(capsys):
    cmd = "bio ncov --fasta --seqid foo --start 10 --end 20"
    run(cmd, capsys=capsys, out="fasta-start.fa")

def test_12(capsys):
    cmd = "bio ncov --fasta --type CDS"
    run(cmd, capsys=capsys, out="CDS.fa")

def test_13(capsys):
    cmd = "bio ncov --fasta --type gene --end 10"
    run(cmd, capsys=capsys, out="gene-start.fa")

def test_14(capsys):
    cmd = "bio ncov --protein --start -10"
    run(cmd, capsys=capsys, out="protein-end.fa")

def test_15(capsys):
    cmd = "bio ncov --translate --type CDS"
    run(cmd, capsys=capsys, out="translate.fa")

def test_16(capsys):
    cmd = "bio ncov:S --fasta --protein --seqid foo"
    run(cmd, capsys=capsys, out="s_prot_foo.fa")

def test_17(capsys):
    cmd = "bio MN996532 --fetch --rename ratg13 --seqid ratg13"
    run(cmd, capsys=capsys, out=None)

def test_18(capsys):
    cmd = "bio ncov ratg13 --end 200 --align"
    run(cmd, capsys=capsys, out="align-dna.txt")

def test_19(capsys):
    cmd = "bio GGC -i --translate"
    run(cmd, capsys=capsys, out="dna-translate.txt")

def test_20(capsys):
    cmd = "bio ncov:S ratg13:S --end 80 --align"
    run(cmd, capsys=capsys, out="align-dna-s.txt")

def test_21(capsys):
    cmd = "bio ncov:S ratg13:S --protein --align"
    run(cmd, capsys=capsys, out="align-protein-s.txt")

def test_22(capsys):
    cmd = "bio ncov:S ratg13:S --end 80 --translate --align"
    run(cmd, capsys=capsys, out="align-translated-s.txt")

def test_23(capsys):
    cmd = "bio THISLINE ISALIGNED  -i --align --local"
    run(cmd, capsys=capsys, out="align-local.txt")

def test_24(capsys):
    cmd = "bio THISLINE ISALIGNED -i --align --global"
    run(cmd, capsys=capsys, out="align-global.txt")

def test_25(capsys):
    cmd = "bio THISLINE ISALIGNED -i --align --semiglobal"
    run(cmd, capsys=capsys, out="align-semiglobal.txt")

def test_26(capsys):
    cmd = "bio 9606 --taxon"
    run(cmd, capsys=capsys, out="taxon_default.txt")

def test_27(capsys):
    cmd = "bio 9606 --lineage --taxon"
    run(cmd, capsys=capsys, out="taxon_lineage.txt")

def test_28(capsys):
    cmd = "bio 9606 --lineage --flat --taxon"
    run(cmd, capsys=capsys, out="taxon_flat_lineage.txt")

def test_29(capsys):
    cmd = "bio ratg13 ncov 9606 --taxon"
    run(cmd, capsys=capsys, out="taxon_mixed.txt")

