# bio: install {#install}

## Prerequisites

   We don't fully automate the dependency installation to give users the option of using the approach they prefer.

To install the requirements with `conda` write:

    conda install -c bioconda biopython pysam python-parasail 
    
You may also install the requirements via `pip`:

    pip install biopython pysam python-parasail

## Installation
    
Install the package with:

    pip install bio --upgrade

## Usage

Type `bio` followed by one or more accession numbers followed by one or more flags or options.

    bio ACC1 [ACC2 ACC3] --option1 value1 --flag1 ...
    
1. A `flag` is a parameter that does not take additional values: `--fetch`
1. An `option` is a parameter that takes an additional value: `--start 100`
    
Use the `-v` flag to produce verbose outputs for each command. 

### Subcommands

Certain words may not be accession numbers as they carry additional meaning and trigger 
alternative actions

    bio align ACC1 [ACC2 ACC3] --option1 value1 --flag1 ...

You may get help on alignments with:

    bio align -h
    
## Getting help

    bio -h
    
to get help for a specific command:
        
        
## Help pages

Besides the default actions `bio` may also take subcommands such as `align`. Each subcommand
has its own command line help page.

#### 1\. Help page for default actions

```{bash, comment=NA}
bio -h
```

#### 2\. Help page for alignments 

```{bash, comment=NA}
bio align -h
```

