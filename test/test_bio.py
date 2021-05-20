
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
    run("bio convert genomes.gb  --fasta > genomes.fa")

def test_4(capsys):
    run("bio convert genomes.gb  --gff > genomes.gff")

