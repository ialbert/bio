
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from biorun.test.generate import *


def test_1(capsys):
    cmd = "bio SARS2 "
    run(cmd, capsys=capsys, out="SARS2.json")


def test_2(capsys):
    cmd = "bio SARS2 --fasta "
    run(cmd, capsys=capsys, out="SARS2.fa")


def test_3(capsys):
    cmd = "bio SARS2 --gff "
    run(cmd, capsys=capsys, out="SARS2.gff")


def test_4(capsys):
    cmd = "bio SARS2 --match ORF1ab --type gene "
    run(cmd, capsys=capsys, out="match.json")


def test_5(capsys):
    cmd = "bio SARS2 --gff --gene S "
    run(cmd, capsys=capsys, out="gene.gff")


def test_6(capsys):
    cmd = "bio SARS2 --gff  --start 10000 --end 20000 "
    run(cmd, capsys=capsys, out="overlap.gff")


def test_7(capsys):
    cmd = "bio SARS2 --gff  --type CDS "
    run(cmd, capsys=capsys, out="type.gff")


def test_8(capsys):
    cmd = "bio SARS2 --fasta --seqid foo --start 10 --end 20 "
    run(cmd, capsys=capsys, out="fasta-start.fa")


def test_9(capsys):
    cmd = "bio SARS2 --fasta --type CDS "
    run(cmd, capsys=capsys, out="CDS.fa")


def test_10(capsys):
    cmd = "bio SARS2 --fasta --type gene --end 10 "
    run(cmd, capsys=capsys, out="gene-start.fa")


def test_11(capsys):
    cmd = "bio SARS2 --protein --start -10 "
    run(cmd, capsys=capsys, out="protein-end.fa")


def test_12(capsys):
    cmd = "bio SARS2:S --fasta --protein --seqid foo "
    run(cmd, capsys=capsys, out="s_prot_foo.fa")


