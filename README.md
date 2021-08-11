# bio: making bioinformatics fun again

> **NOTE**: Exploratory work in progress. Major changes may be introduced at any time.

`bio` - command-line utilities to make bioinformatics explorations more enjoyable.

Typical usage:

    bio task --param data

where task can be: `fetch`, `fasta`, `align` and many others.

## Quick links

* Documentation: https://www.bioinfo.help/
* Usage examples: [usage.sh][usage]

[docs]: https://ialbert.github.io/bio/
[usage]: https://github.com/ialbert/bio/blob/master/test/usage.sh

## What does this software do?

This software was designed to teach bioinformatics concepts.

If you've ever done bioinformatics, you know how even seemingly straightforward tasks require multiple steps, arcane incantations, and various other preparations that slow down progress. 

Even well-defined, supposedly simple tasks can take a seemingly inordinate number of complicated steps. The `bio` package is meant to solve that tedium. Suppose we fetch multiple accession numbers from NCBI into a genbank file.

    bio fetch NC_045512,MN996532 > genomes.gb

Now that we have GenBank file with multiple genomes, but how do we extract the information from them? This is where we use `bio`

    bio fasta genomes.gb --end 10

prints:

    >NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome [1:10]
    ATTAAAGGTT
    >MN996532.2 Bat coronavirus RaTG13, complete genome [1:10]
    ATTAAAGGTT

whereas to get the coding sequence coordinate that correspond to gene `S` we can write:
 
    bio gff genomes.gb --gene S 

prints:

    ##gff-version 3
    NC_045512.2     .       CDS     21562   25384   .       +       .       ID=1;Name=YP_009724390.1;Parent=YP_009724390.1
    MN996532.2      .       CDS     21559   25369   .       +       .       ID=2;Name=QHR63300.2;Parent=QHR63300.2

## `bio` is stream oriented

`bio` supports stream oriented programming where the output of one task may be chained into the second. Take the example above
but now start with a file `acc.txt` that contains just the accession numbers:

    NC_045512
    MN996532

we can run `bio` to find the first three codons for each coding sequence for gene `S`:

     cat acc.txt | bio fetch | bio fasta --gene S --end 9

to print:

    >YP_009724390.1 CDS surface glycoprotein [1:9]
    ATGTTTGTT
    >QHR63300.2 CDS spike glycoprotein [1:9]
    ATGTTTGTT

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

Detailed documentation is maintained at

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
