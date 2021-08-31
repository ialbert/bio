
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from generate import *

# Initialize directories.
init_dirs()


def test_1(capsys):
    run("set -uex")

def test_2(capsys):
    run("bio fetch NC_045512 MN996532 > genomes.gb")

def test_3(capsys):
    run("echo NC_045512 | bio fetch > sars2.gb")

def test_4(capsys):
    run("bio fasta genomes.gb --end  100 > genomes.fa")

def test_5(capsys):
    run("cat genomes.gb | bio fasta --end  100 > genomes.fa")

def test_6(capsys):
    run("bio fasta genomes.gb --end  100  --alias alias.txt > genomes.alias.fa")

def test_7(capsys):
    run("bio fasta genomes.gb --end 10 --type CDS > cds.fa")

def test_8(capsys):
    run("bio fasta genomes.gb --type CDS --translate > translate.fa")

def test_9(capsys):
    run("bio fasta genomes.gb  --protein > protein.fa")

def test_10(capsys):
    run("bio fasta -s -3 > stop.fa")

def test_11(capsys):
    run("bio align GATTACA GATCA > gattaca1.txt")

def test_12(capsys):
    run("bio align GATTACA GATCA --global > gattaca2.txt")

def test_13(capsys):
    run("bio align GATTACA GATCA --local > gattaca3.txt")

def test_14(capsys):
    run("bio align GATTACA GATCA --variant > variant.txt")

def test_15(capsys):
    run("bio fasta --gene S --protein  genomes.gb > s.fa")

def test_16(capsys):
    run("bio align s.fa > align-default.txt")

def test_17(capsys):
    run("bio align s.fa --table > align-table.txt")

def test_18(capsys):
    run("bio align s.fa --variant > align-variant.txt")

def test_19(capsys):
    run("bio gff genomes.gb > genomes.gff")

def test_20(capsys):
    run("bio gff genomes.gb --type CDS > CDS.gff")

def test_21(capsys):
    run("bio gff -s 300 -e 10k > slice.gff")

def test_22(capsys):
    run("bio taxon 117565 -d 5 > taxonomy.txt")

def test_23(capsys):
    run("bio taxon genomes.gb --lineage > lineage.txt")

def test_24(capsys):
    run("bio data 11138 -H > meta.txt")

def test_25(capsys):
    run("bio explain exon > so.txt")

def test_26(capsys):
    run("bio explain food vacuole > go.txt")

def test_27(capsys):
    run("bio explain neutral > search.txt")

