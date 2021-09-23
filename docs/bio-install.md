# How to install bio {#bio-install}

`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. 

## Installation

Installation:

    pip install bio --upgrade

## Download prebuild databases

Some tasks (`bio taxon`, `bio explain`) need access to additional data. Install the required databases with:

    bio --download

## Test the code

To execute a batch of tests run

    bio test

## Usage

Type `bio` followed by a task (`fasta`, `gff`, `taxon`) then followed by one or more flags or options.

    bio fetch NC_045512 > genome.gb
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
