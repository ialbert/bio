# bio: making bioinformatics fun again

> The software is currently under development. It is operational but not fully vetted.

`bio` - command-line utilities to make bioinformatics explorations more enjoyable.

## Quick links

* Documentation: https://ialbert.github.io/bio/
* Usage examples: [bio_examples.sh][examples]

[docs]: https://ialbert.github.io/bio/
[examples]: https://github.com/ialbert/bio/blob/master/test/bio_examples.sh

## Why do we need this software?

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, reading documentation and numerous other preparations that slow down your progress. 

Time and again I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package is meant to solve that tedium.  With `bio` you can write things like this:

    # Fetch the data from NCBI.
    bio NC_045512 --fetch --rename ncov
    bio MN996532  --fetch --rename ratg13
    
    # Align the DNA for the S protein.
    bio ncov:S ratg13:S --end 90 --align

to align the first 90 basepairs of the DNA sequence of the  `S` protein from the SARS-COV-2 novel coronavirus to its closest (known) relative, the bat coronavirus RaTG13. The command above will print:

```
### 1: YP_009724390 vs QHR63300.2 ###

Length: 90 (semiglobal)
Query:  90 [1, 90]
Target: 90 [1, 90]
Score:  387
Ident:  83/90 (92.2%)
Simil:  83/90 (92.2%)
Gaps:   0/90 (0.0%)
Matrix: nuc44(-11, -1)

YP_009724390 ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAAT
           1 ||||||||||||||||||||||||||||||||.||||||||||||||||||||.|||||.||||||||.|||||.|||||||||||.||. 90
QHR63300.2   ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTTTCTAGTCAGTGTGTTAATCTAACAACTAGAACTCAGTTACCTCCTGCATACACCAAC
```

If you wanted to align the same sequences as translated proteins `bio` lets you write:

    bio ncov:S ratg13:S --end 90 --translate --align
    
to generate:

```
### 1: YP_009724390 vs QHR63300.2 ###

Length: 30 (semiglobal)
Query:  30 [1, 30]
Target: 30 [1, 30]
Score:  153
Ident:  30/30 (100.0%)
Simil:  30/30 (100.0%)
Gaps:   0/30 (0.0%)
Matrix: blosum62(-11, -1)

YP_009724390 MFVFLVLLPLVSSQCVNLTTRTQLPPAYTN
           1 |||||||||||||||||||||||||||||| 30
QHR63300.2   MFVFLVLLPLVSSQCVNLTTRTQLPPAYTN
```

Beyond alignments there is a lot more to `bio` and we recommend looking at the [documentation][docs]

## Who is `bio` designed for?

The software was written to teach bioinformatics and is the companion software to the [Biostar Handbook][handbook] textbook. The targeted audience comprises:

- Students learning about bioinformatics.
- Bioinformatics educators that need a platform to demonstrate bioinformatics concepts. 
- Scientists working with large numbers of similar genomes (bacterial/viral strains).
- Scientists that need to closely investigate and understand particular details of a genomic region.

The ideas and motivations fueling the creation of `bio` came to us while educating the many cohorts of students that used the handbook in the classrom. 

You see, in bioinformatics, many tasks that should be straightforward are, instead, needlessly complicated. `bio` is an opinionated take on how bioinformatics, particularly data presentation and access, should be simplified. 

[handbook]: https://www.biostarhandbook.com/

## Documentation

The documentation is maintained at

* https://ialbert.github.io/bio/


## Quick install
    
`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. While the following command should work:

    pip install bio --upgrade
        
we usually recommend installing prerequisites with conda like so:

    conda install -c bioconda biopython parasail-python
    
then proceed with:

    pip install bio --upgrade

    
See more details in the [documentation][docs].

## Development

If you clone the repository we recommend to install as development package with:

    python setup.py develop
    
## Testing

Testing uses the pytest framework:

    pip install pytest

To run all tests use:

    make test
    
Tests are automatically built from a test script that mimics real life usage scenarios.

* https://github.com/ialbert/bio/blob/master/test/bio_examples.sh

## New tests

To add a new test first run the command you wish to test, for example:

    bio foo --gff > output.gff

in the `test/data` directory. After that add the same command above into the master script:

* https://github.com/ialbert/bio/blob/master/test/bio_examples.sh
    
followed by:

    make build_tests
    
The latter command will automatically generate a Python test for each line in the master script.

The automatically generated test will verify that the command is operational and that the output matches the expectations.
