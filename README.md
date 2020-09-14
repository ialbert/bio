# Bioinformatics Utilities

Making bioinformatics great again.

Series of utility functions that streamline bioinformatics.


## How to use

> **Note**: data with NCBI accession numbers will be automatically downloaded from the corresponding data repositories then stored locally. No need to obtain data beforehand.

Align the DNA sequence of the `S` protein of two two `SARS-COV-2` genomes:

    bio align ACC1:S ACC2:S

It prints:

    ------

Align the `S` protein as peptide sequences:

    bio align --protein ACC1:S ACC2:S


## Requirements

The package requires `BioPython` and `pysam`. Depending on your work style use:

    conda install biopython pysam
    
or

    pip install biopython pysam
    
    
## Installation
        
    pip install bio
    
    
## Development mode

    git clone https://github.com/ialbert/bio.git
    (cd bio && python setup.py develop)
    

## Run

Type

    bio 
        
to see the usage:

    Bioinformatics utilities: 0.0.1

    Usage: bio COMMAND
    
    Data commands:
    
        convert    - convert biological data to other formats
    
    Operations:
    
        align      - align sequences with different algorithms
    
    Get more help on each command with:
    
        bio COMMAND -h    
