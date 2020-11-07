# bio: install {#install}

## Prerequisites

To install the requirements write:

    conda install biopython pysam python-parasail 
    
You may also install the requirements via `pip`:

    pip install biopython pysam python-parasail
   
We don't fully automate the dependency installation to give users the option of using the approach they prefer.

## Installation
    
Install the package with:

    pip install bio --upgrade

## Usage

Type `bio` followed by one or more accession numbers followed by one or more flags or options.

    bio ACC1 [ACC2 ACC3] --option1 --option2 --flag1 ...
    
1. A `flag` is a parameter that does not take additional values: `--fetch`
1. An `option` is a paramter that takes an additional value: `--start 100`
    
Use the `-v` command to produce verbose outputs for each command. 

## Getting help

    bio -h
    
to get help for a specific command:
        
## Example help page

```{bash, comment=NA}
bio -h
```


