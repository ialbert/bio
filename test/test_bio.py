
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
    run("bio fetch NC_045512 MN996532 --quiet > genomes.gb")

def test_3(capsys):
    run("bio convert genomes.gb  --fasta > genomes.fa")

def test_4(capsys):
    run("bio convert genomes.gb --end  10 > slice.fa")

def test_5(capsys):
    run("bio convert genomes.gb --end 10 --features > features.fa")

def test_6(capsys):
    run("bio convert genomes.gb --end 10 --type CDS > cds.fa")

def test_7(capsys):
    run("bio convert genomes.gb --type CDS --translate > translate.fa")

def test_8(capsys):
    run("bio convert genomes.gb  --protein > protein.fa")

def test_9(capsys):
    run("bio convert genomes.gb  --gff > genomes.gff")

