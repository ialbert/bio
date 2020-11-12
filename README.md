# Introduction

> Under development (the package is functional but not fully vetted)

**Making bioinformatics fun again.**

`bio` - command-line utilities to make bioinformatics explorations more enjoyable.

## Why do we need this software?

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, reading documentation and numerous other preparations that slow down your progress. 

Time and again I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package is meant to solve that tedium.  With `bio` you can write this do this:

    # Fetch the data
    bio NC_045512 --fetch --rename ncov
    bio MN996532  --fetch --rename ratg13
    
    # Run an alignment.
    bio align ncov:S ratg13:S --end 80 --translate

to align the first 80 basepairs DNA sequence translated as for the proteins taken from SARS-COV-2 and its closest (known) relative Bat coronavirus RatG13. Read more about the process in the documentation.

## Learn more about how `bio` works

The documentation is maintained at

* https://ialbert.github.io/bio/

Or in the github repository as markdown files:

* https://github.com/ialbert/bio/tree/master/md

## Quick install
    
Install the package with:

    # Prerequisites
    conda install -c bioconda biopython pysam parasail-python
    
    # Install the bio package.
    pip install bio --upgrade
    
See more details in the documentation.

## Development

If you clone the repository we recommend to install as development package with:

    python setup.py develop
    
## Testing

To run all tests use:

    make test
    
Tests are automatically built from a test script that mimics real life usage scenarios.

* https://github.com/ialbert/bio/blob/master/test/test_bio_data.sh

To add a new test first run the command you wish to test

    bio foo --gff > output.gff

in the `test/data` directory. After that add the same command above into the master script:

* https://github.com/ialbert/bio/blob/master/test/test_bio_data.sh
    
followed by:

    make build_tests
    
The latter command will automatically generate a Python test for each line in the master script.
The automatically generated test will verify that the command is operational and that the output matches the expectations.
