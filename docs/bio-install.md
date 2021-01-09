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

    bio fetch NC_045512 -v

now list the known data:

    bio data

try out a conversion:

    bio convert NC_045512

## Usage

Type `bio` followed by one or more accession numbers or data names followed by one or more flags or options.

## Getting help

    bio -h
    
## Tasks

Tasks will have separate help pages that are shown when the command is run with no parameters (or using the `-h` flag)
