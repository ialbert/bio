# How to install bio {#bio-install}

`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem.

## Contact information

If you are reporting a bug please open an issue below:

*  https://github.com/ialbert/bio/issues

If you have other questions/suggestions on how to use  `bio` use the discussion board:

* https://github.com/ialbert/bio/discussions

## Installation

Installation:

    pipx install bio --upgrade

or with:

    pip install bio --upgrade

Some tasks (`bio taxon`, `bio explain`) need access to additional data. Install the required databases with:

    bio --download

To execute a batch of tests run

    bio test

Tests depend on the availability of a number of online services. Some test fail whenever the corresponding remote data service is unavailable. Connection problems are typically corrected within a short period of time.

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

Use the discussion board for other questions:

* https://github.com/ialbert/bio/discussions

## Tasks

Tasks will have separate help pages that are shown when the command is run with no parameters (or using the `-h` flag)
