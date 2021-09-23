# bio: general tips {#bio-tips}

Substantial effort has been devoted to making the command line more error-tolerant and user friendly.

## See the help

Each command, when invoked with no parameters, will produce help:

    bio

## All commands can read the standard input

    bio fetch AF086833  KJ660346 > genome.gb

or if the file `ids.txt` contains:

    AF086833
    KJ660346

you can run:

    cat ids.txt | bio fetch > genome.gb

## Some commands can create sequences on the fly

    bio align GATCA GATTACA

or

    bio fasta GATTACA --translate

## Parameter order does not matter

You may write:

    bio fasta --type CDS --end 10 genomes.gb

or:

    cat genomes.gb | bio fasta --type CDS --end 10

both will work and produce the same results.

## Parameter action

Each parameter will be applied sequentially in an internally determined order that makes the most sense: 

    bio fasta --type CDS --end 10  --translate genomes.gb

will produce the same results as:

    bio fasta --translate --end 10 --type CDS  genomes.gb
    
Both commands first select `CDS` types, apply a slice on each sequence, and then use the translation operator.
 
## Multiple accession numbers
   
Many commands allow using multiple accession numbers; in that case, the operations will take place sequentially on each.

    bio fetch NC_045512 MN996532

or comma separated:

    bio fetch NC_045512,MN996532
  
## Parameter forms

You may use single or double dashes on parameters:

    bio fasta genomes.gb --end 100
    
The command above is equivalent to:

    bio fasta genomes.gb -end 100


## The coordinate system is 1 based

Coordinates are one based (inclusive on both ends) identical to GFF coordinate formats.

    bio fasta --start 10 --end 20 genomes.gb
    
The interval of 10 to 20 is 11 bases long! To make a single base long slice start and end on the same value:

    bio fasta --start 10 --end 10 genomes.gb

## Number formatting

Numbers for start and end coordinates may be written in human-friendly forms, like so: 

    bio fasta -start 1kb -end 2kb genomes.gb

accepted formats:

* `5000` 
* `5,000`
* `5k` or `5kb`
* `5K` or `5KB`  

