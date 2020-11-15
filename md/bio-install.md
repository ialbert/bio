# Software installation {#install}

## Prerequisites

We don't fully automate the dependency installation to give users the option of using the approach they prefer. To install the requirements with `conda` write:

    conda install -c bioconda biopython pysam parasail-python
    
You may also install the requirements via `pip`:

    pip install biopython pysam parasail-python

## Installation
    
Install the package with:

    pip install bio --upgrade

## Quick start

Run a simple fetch command, set the verbose mode to see what is happening:

    bio NC_045512 --fetch -v
 
now list the known data:

    bio --list
    
try out a conversion:

    bio NC_045512 --gff
    
## Usage

Type `bio` followed by one or more accession numbers followed by one or more flags or options.

    bio [command] [words] --option value --flag 
    
1. Commands may be: `align`, `taxon`.  When no commands are passed the default actions take place.
1. The `words` may be one or more data names.    
1. A `flag` is a parameter that does not take additional values: `--fetch` or `--list`
1. An `option` is a parameter that takes an additional value: `--start 100` or `--gap-open 10`
    
Use the `-v` flag to produce verbose outputs for a run. 

### Subcommands

Certain words may not be used as data names  as they carry additional meaning and trigger 
alternative actions

    bio align data1 data2 --option value --flag

You may get help on alignments with:

    bio align -h
    
## Getting help

    bio -h
    
to get help for a specific command:
        
        
## Help pages

Besides the default actions `bio` may also take subcommands such as `align`. Each subcommand
has its own command line help page.

### 1\. Help page for default actions

```{bash, comment=NA}
bio -h
```

### 2\. Help page for alignments 

```{bash, comment=NA}
bio align -h
```

### 3\. Help page for taxonomy 

```{bash, comment=NA}
bio taxon -h
```

