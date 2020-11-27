# Software installation {#install}

`bio` works on Linux and Mac computers and on Windows when using the Linux Subsystem. 

## Installation

While the following command should work:

    pip install bio --upgrade

we usually recommend installing prerequisites with conda like so:

    conda install -c bioconda biopython parasail-python

then proceed with:

    pip install bio --upgrade

## Quick start

Run a simple fetch command, set the verbose mode to see what is happening:

    bio NC_045512 --fetch -v

now list the known data:

    bio --list

try out a conversion:

    bio NC_045512 --gff

## Usage

Type `bio` followed by one or more accession numbers or data names followed by one or more flags or options.

## Subcommands

Certain flags trigger different behaviors:

    bio --align
    bio --taxon
    bio --sra 
    
## Getting help

    bio -h

Subcommands will have separate help pages. For example:

    bio --align -h

## Actual help pages

Below we include the help page for each active command

### 1\. Help page for default actions

```{bash, comment=NA}
bio -h
```

### 2\. Help page for alignments 

```{bash, comment=NA}
bio --align -h
```

### 3\. Help page for taxonomy 

```{bash, comment=NA}
bio --taxon -h
```

### 4\. Help page for SRA search 

```{bash, comment=NA}
bio --sra -h
```
