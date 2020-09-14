# bio: fetch {#fetch}

Obtains data identifed via an accession number from NCBI (or EMBL):

`fetch` will store the downloaded data in a cache directory `~/bio` so that next time around it does not need to connect to the internet.

All other commands (`convert`, `align` etc) that can operate via accession numbers also use the `fetch` behind the scences.

## Examples

    # Get a single accession number
    bio fetch NC_045512 > results.gb

    # Get multiple accession numbers into a single file.
    bio fetch NC_045512 AF086833 > results.gb

    # Extract the S protein's DNA sequence as FASTA
    bio fetch NC_045512 --name S > nucleotide.fa

    # Extract the S protein's DNA sequence as FASTA
    bio fetch NC_045512 --name S --range 1-200 > nucleotide.fa
    
    # Extract the S protein's translated sequence from the GenBank file
    bio fetch NC_045512:S --trans > protein.fa

    # Extract the S protein's DNA coordinates from the GenBank file
    bio fetch NC_045512:S --gff > results.gff

## Usage
   
```{bash, comment=NA}
bio fetch NC_045512 | head
```

## Command line help
```{bash, comment=NA}
bio fetch -h
```
