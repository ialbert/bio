# bio: making bioinformatics fun again

> **NOTE**: Exploratory work in progress. Major changes may be introduced at any time.

`bio` - command-line utilities to make bioinformatics explorations more enjoyable.

Typical usage:

    bio task data --parameter

where task can be: `fetch`, `convert`, `align` and many others.

## Quick links

* Documentation: https://www.bioinfo.help/
* Usage examples: [usage.sh][usage]

[docs]: https://ialbert.github.io/bio/
[usage]: https://github.com/ialbert/bio/blob/master/test/usage.sh

## What does this software do?

This software was designed to teach bioinformatics concepts.

If you've ever done bioinformatics, you know how even seemingly straightforward tasks require multiple steps, arcane incantations, and various other preparations that slow down progress. 

Even well-defined, supposedly simple tasks can take a seemingly inordinate number of complicated steps. The `bio` package is meant to solve that tedium.  With `bio`, you can write things like this:

    # Fetch multiple accession numbers from genbank
    bio fetch NC_045512 MN996532 > genomes.gb

    # Convert Genbank to FASTA.
    bio convert genomes.gb  --fasta

    # Convert Genbank to GFF.
    bio convert genomes.gb  --gff
    

## Who is `bio` designed for?

The software was written to teach bioinformatics and is the companion software to the [Biostar Handbook][handbook] textbook. The targeted audience comprises:

- Students learning about bioinformatics.
- Bioinformatics educators who need a platform to demonstrate bioinformatics concepts. 
- Scientists working with large numbers of similar genomes (bacterial/viral strains).
- Scientists who need to investigate and understand the precise details of a genomic region closely.

The ideas and motivations fueling `bio` have been developed while educating the many cohorts of students who used the handbook in the classroom. 

They have all noted that when practicing bioinformatics, many tasks that should be straightforward are, instead, needlessly complicated. `bio` is an opinionated take on how bioinformatics, particularly data representation and access, should be simplified and streamlined.

[handbook]: https://www.biostarhandbook.com/

## Documentation

The documentation is maintained at

* https://www.bioinfo.help/

## Quick install
    
`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. 

    pip install bio --upgrade
            
See more details in the [documentation][docs].

## Development

If you clone the repository, we recommend that you install it as a development package with:

    python setup.py develop
    
## Testing

Testing uses the `pytest` framework:

    pip install pytest

To run all tests, use:

    make test
    
Tests are automatically built from a shell script that mimics real-life usage scenarios.

* https://github.com/ialbert/bio/blob/master/test/usage.sh


## Generating documentation

To generate the docs, you will need the `bookdown` package:

    conda install r-bookdown r-servr
    
To run the docs in a browse:
    
    make 
    
then visit http://localhost:8000

To render the docs write:

    make docs

To push out the latest docs:    
    
    make sync
    
