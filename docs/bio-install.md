# Installation {#bio-install}

`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. 

## Installation

While the following command should work:

    pip install bio --upgrade

We usually recommend installing prerequisites with conda like so:

    conda install -c bioconda biopython

then proceed with `pip install bio --upgrade`.

## Quick start

Run a simple fetch command, set the verbose mode to see what is happening:

    bio NC_045512 --fetch -v

now list the known data:

    bio --list

try out a conversion:

    bio NC_045512 --gff

## Usage

Type `bio` followed by one or more accession numbers or data names followed by one or more flags or options.

## Getting help

    bio -h
    
## Subcommands

Certain flags trigger different functionality. We call these flags subcommands.

Subcommands will have separate help pages that are shown when the command is run with no parameters (or using the `-h` flag)

    bio --fasta
    bio --align
    bio --taxon
    bio --defint
    bio --sra
