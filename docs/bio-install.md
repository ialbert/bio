# Installation {#bio-install}

`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. 

## Installation

We usually recommend installing prerequisites with conda like so:

    conda install -c bioconda biopython

Then following up with:

    pip install bio --upgrade

## Usage

Type `bio` followed by a task (`fasta`, `gff`, `taxon`) then followed by one or more flags or options.

    bio fasta genome.gb
    bio gff genome.gb --type CDS
    bio taxon 2697049 --lineage

## Getting help

You can get help on any command by just entering it

    bio 
    bio fasta
    bio taxon

each line above will print a help on the task usage.

## Tasks

Tasks will have separate help pages that are shown when the command is run with no parameters (or using the `-h` flag)
