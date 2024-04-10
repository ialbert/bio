# bio: making bioinformatics fun again

`bio` - command-line utilities to make bioinformatics explorations more enjoyable.

`bio` is a bioinformatics toy to play with.

Like LEGO pieces that match one another `bio` aims to provide you with commands that naturally fit together and let you express your intent with short, explicit and simple commands. It is a project in an exploratory phase, we'd welcome input and suggestions on what it should grow up into.

## What does this software do?


If you've ever done bioinformatics, you know how even seemingly straightforward tasks require multiple steps, arcane incantations, and various other preparations that slow down progress. 

Even well-defined, supposedly simple tasks can take a seemingly inordinate number of complicated steps. The `bio` package is meant to solve that tedium. 

## Usage examples

    # Fetch genbank data
    bio fetch NC_045512 MN996532 > genomes.gb

    # Convert the first then bases of the genomes to FASTA.
    bio fasta genomes.gb --end 10

    # Align the coding sequences for the S protein
    bio fasta genomes.gb --gene S --protein | bio align | head

    # Print the GFF record that corresponds to the coding sequence for gene S
    bio gff genomes.gb --gene S 

    # Show the descendants of taxid 117565
    bio taxon 117565 | head

    # Show the lineage of a taxonomic rank.
    bio taxon 117565 --lineage | head

    # Get metadata on a viral sample
    bio meta 11138 -H | head

    # Define a sequence ontology terms
    bio define exon

    # Define a gene ontology terms
    bio define food vacuole

## Documentation

Detailed documentation is maintained at

* https://www.bioinfo.help/

## Quick install
    
`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. 

As a rule, all Python based command line utilities should be installed via [pipx][pipx] to avoid conflicts with other Python packages:

[pipx]: https://pipx.pypa.io/stable/

    pipx install bio 

Alternatively, if you can also use `pip` to install:

    pip install bio 
            
See more details in the [documentation][docs].

## `bio` is stream oriented

`bio` supports stream oriented programming where the output of one task may be chained into the second. Take the example above
but now start with a file `acc.txt` that contains just the accession numbers:

    NC_045512
    MN996532

we can run `bio` to generate a VCF file with the variants of the S nucleotides forming the S protein like so:

    cat acc.txt | bio fetch | bio fasta --gene S | bio align --vcf | head

to print:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">
    ##contig=<ID=YP_009724390.1,length=3822,assembly=YP_009724390.1>
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  QHR63300.2
    YP_009724390.1  33      33C/T   C       T       .       PASS    TYPE=SNP        GT      1
    YP_009724390.1  54      54T/A   T       A       .       PASS    TYPE=SNP        GT      1
    YP_009724390.1  60      60C/T   C       T       .       PASS    TYPE=SNP        GT      1
    YP_009724390.1  69      69A/G   A       G       .       PASS    TYPE=SNP        GT      1


## Who is `bio` designed for?

The software was written to teach bioinformatics and is the companion software to the [Biostar Handbook][handbook] textbook. The targeted audience comprises:

- Students learning about bioinformatics.
- Bioinformatics educators who need a platform to demonstrate bioinformatics concepts. 
- Scientists working with large numbers of similar genomes (bacterial/viral strains).
- Scientists who need to investigate and understand the precise details of a genomic region closely.

The ideas and motivations fueling `bio` have been developed while educating the many cohorts of students who used the handbook in the classroom. `bio` is an opinionated take on how bioinformatics, particularly data representation and access, should be simplified and streamlined.

[handbook]: https://www.biostarhandbook.com/
[docs]: https://www.bioinfo.help/

## Development

We use the `hatch` build system to manage the software:

https://github.com/pypa/hatch

You can either use `hatch` or `pip` to install the software in editable mode:

    pip install --editable .
    
## Testing

`bio` can test itself, to run all tests execute:

    bio test

Tests are automatically built from a shell script that mimics real-life usage scenarios.

* https://github.com/ialbert/bio/blob/master/src/biorun/data/usage.sh

## Generating documentation

To generate the docs, you will need the `bookdown` package:

    conda install r-bookdown r-servr
    
To run the docs in a browse:
    
    make 
    
then visit http://localhost:8000
