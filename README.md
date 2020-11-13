# bio: making bioinformatics fun again

> Under development (the package is functional but not fully vetted)

`bio` - command-line utilities to make bioinformatics explorations more enjoyable.

## Why do we need this software?

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, reading documentation and numerous other preparations that slow down your progress. 

Time and again I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package is meant to solve that tedium.  With `bio` you can write things like this:

    # Fetch the data from NCBI.
    bio NC_045512 --fetch --rename ncov
    bio MN996532  --fetch --rename ratg13
    
    # Align the DNA for the S protein.
    bio align ncov:S ratg13:S --end 90 

to align the first 90 basepairs of the DNA sequence of the `S` protein,  taken from SARS-COV-2 and its closest (known) relative bat coronavirus RatG13. If you wanted to align the sequences as translated proteins you would write:

    bio align ncov:S ratg13:S --end 90 --translate
    
but just to make sure, there is a lot more to `bio` than alignments.

## Who is `bio` designed for?

- Students learning about bioinformatics.
- Bioinformatics educators that need a platform to demonstrate bioinformatics concepts. 
- Scientists working with large numbers of similar genomes (bacterial/viral strains).
- Scientists that need to closely investigate and understand particular details of a genomic region.

`bio` is designed for detail oriented investigations. 

## Learn more about how `bio` works

The documentation is maintained at

* https://ialbert.github.io/bio/

## Quick install
    
Install the package with:

    # Install prerequisites with conda, pip should also work.
    conda install -c bioconda biopython pysam parasail-python
    
    # Install the bio package.
    pip install bio --upgrade
    
See more details in the documentation.

## Development

If you clone the repository we recommend to install as development package with:

    python setup.py develop
    
## Testing

Testing uses the pytest framework:

    pip install pytest

To run all tests use:

    make test
    
Tests are automatically built from a test script that mimics real life usage scenarios.

* https://github.com/ialbert/bio/blob/master/test/test_bio_data.sh

## Adding new tests

To add a new test first run the command you wish to test

    bio foo --gff > output.gff

in the `test/data` directory. After that add the same command above into the master script:

* https://github.com/ialbert/bio/blob/master/test/test_bio_data.sh
    
followed by:

    make build_tests
    
The latter command will automatically generate a Python test for each line in the master script.
The automatically generated test will verify that the command is operational and that the output matches the expectations.
