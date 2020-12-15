# bio: making bioinformatics fun again

> The software is currently under development. It is operational but not well tested and evaluated.

`bio` - command-line utilities to make bioinformatics explorations more enjoyable.

## Quick links

* Documentation: https://www.bioinfo.help/
* Usage examples: [bio-examples.sh][examples]

[docs]: https://ialbert.github.io/bio/
[examples]: https://github.com/ialbert/bio/blob/master/test/bio-examples.sh

## What does this software do?

This software was designed to teach bioinformatics concepts.

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, and various other preparations that slow down progress. 

Even well defined, supposedly simple tasks can take seemingly inordinate number of complicated steps. The `bio` package is meant to solve that tedium.  With `bio` you can write things like this:

    bio NC_045512 --fetch --rename ncov
    bio MN996532  --fetch --rename ratg13
    
to fetch the data from NCBI and rename data to more meaningful labels, then write:

    bio ncov:S ratg13:S --end 60 --align -pep1

to align the DNA for the S protein while also showing the translation with 1 letter peptide code:

```
# Ident=57(95.0%)  Mis=3(5.0%)  Gaps=0(0.0%)  Target=(1, 60)  Query=(1, 60)  Length=60  Score=273.0  NUC.4.4(11,1)

              M  F  V  F  L  V  L  L  P  L  V  S  S  Q  C  V  N  L  T  T
YP_009724390 ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
             ||||||||||||||||||||||||||||||||.||||||||||||||||||||.|||||. 60
QHR63300.2   ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTTTCTAGTCAGTGTGTTAATCTAACAACT
              M  F  V  F  L  V  L  L  P  L  V  S  S  Q  C  V  N  L  T  T
```

`bio` was designed to use words that make sense: align, translate, complement, protein, taxon, type etc. If you wanted to align the same sequences when translated into proteins `bio` lets you write:

    bio ncov:S ratg13:S --end 60 --translate --align
    
to generate:

```
# Ident=20(100.0%)  Mis=0(0.0%)  Gaps=0(0.0%)  Target=(1, 20)  Query=(1, 20)  Length=20  Score=98.0  BLOSUM62(11,1)

YP_009724390 MFVFLVLLPLVSSQCVNLTT
             |||||||||||||||||||| 20
QHR63300.2   MFVFLVLLPLVSSQCVNLTT
```

Beyond alignments there is a lot more to `bio` and we recommend looking at the [documentation][docs]

## Who is `bio` designed for?

The software was written to teach bioinformatics and is the companion software to the [Biostar Handbook][handbook] textbook. The targeted audience comprises:

- Students learning about bioinformatics.
- Bioinformatics educators that need a platform to demonstrate bioinformatics concepts. 
- Scientists working with large numbers of similar genomes (bacterial/viral strains).
- Scientists that need to closely investigate and understand particular details of a genomic region.

The ideas and motivations fueling the creation of `bio` have been developed while educating the many cohorts of students that used the handbook in the classroom. 

You see, in bioinformatics, many tasks that should be straightforward are, instead, needlessly complicated. `bio` is an opinionated take on how bioinformatics, particularly data representation and access, should be simplified. 

[handbook]: https://www.biostarhandbook.com/

## Documentation

The documentation is maintained at

* https://www.bioinfo.help/

## Quick install
    
`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. 

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

* https://github.com/ialbert/bio/blob/master/test/bio-examples.sh

## New tests

To add a new test first run the command you wish to test, for example:

    bio foo --gff > output.gff

in the `test/data` directory. After that add the same command above into the master script:

* https://github.com/ialbert/bio/blob/master/test/bio-examples.sh
    
followed by:

    make build_tests
    
The latter command will automatically generate a Python test for each line in the master script.

The automatically generated test will verify that the command is operational and that the output matches the expectations.


## Generating documentation

To generate the docs you will need `bookdown`

    conda install r-bookdown r-servr
    
To run the docs in a browse:
    
    make 
    
then visit http://localhost:8000

To render the docs write:

    make docs

To push out the latest docs:    
    
    make sync
    