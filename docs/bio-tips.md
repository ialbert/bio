# General tips {#bio-tips}

Substantial effort has been devoted to making the command line more error-tolerant and user friendly.

## See the help

Each command, when invoked with no parameters, will produce help:

    bio --align 
    
## Print more information

Use the `-v` flag to produce verbose outputs for each command.

## Parameter order 

You may write:

    bio convert ncov --fasta --type CDS --end 10 

or:

    bio convert --type CDS --end 10 --fasta ncov 

both will work and produce the same results.

## Parameter action

Each parameter will be applied sequentially in an internally determined order that makes the most sense: 

    bio convert ncov --fasta --type CDS --end 10  --translate

will produce the same results as:

    bio convert ncov --fasta --translate --end 10 --type CDS 
    
Both commands first select `CDS` types, apply a slice on each sequence, and then use the translation operator.
 
## Multiple accession numbers
   
Many commands allow using multiple accession numbers; in that case, the operations will take place sequentially on each.

    bio fetch NC_045512 MN996532 
  
## Parameter forms

You may use single or double dashes on parameters:

    bio ncov --fasta --end 100
    
The command above is equivalent to:

    bio convert ncov -fasta -end 100
    
## Interactive mode

Passing the `-i` flag allows data to be passed from the command line. For example:

```{bash, comment=NA}
bio convert  ATGATTATATATA --translate -i 
```

Note how the input was read as parameters from the command line. We make use of this feature when explicitly exploring simple data.

## The coordinate system is 1 based

Coordinates are one based (inclusive on both ends) identical to GFF coordinate formats.

    bio convert ncov -fasta -start 10 --end 20
    
The interval of 10 to 20 is 11 bases long! To make a single base long slice start and end on the same value:

    bio convert ncov -fasta -start 10 --end 10

## Number formatting

Numbers for start and end coordinates may be written in human-friendly forms, like so: 

    bio convert ncov -fasta -start 1kb --end 2kb 

accepted formats:

* `5000` 
* `5,000`
* `5k` or `5kb`
* `5K` or `5KB`  

