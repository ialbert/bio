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
    run("bio fetch NC_045512 KM233118 > genome.gb")
