# bio: view

The `view` command can seamlessly:

- convert data from GenBank/EMBL to other formats: FASTA, GFF
- extract sections of data: features by name, coordinate, range
- download data identified via an accession number from NCBI or EMBL.

You can also use `fetch` to

`fetch` will store the downloaded data in a cache directory `~/bio` so that next time around it does not need to connect to the internet.

All other commands (`convert`, `align` etc) that can operate via accession numbers also use the `fetch` behind the scences.

## Examples

    # Get a single accession number
    bio view NC_045512 > results.gb

    # Get multiple accession numbers into a single file.
    bio view NC_045512 AF086833 > results.gb

    # Extract the S protein's DNA sequence as FASTA
    bio view NC_045512 --name S > nucleotide.fa

    # Extract the S protein's DNA sequence as FASTA
    bio view NC_045512 --name S --range 1-200 > nucleotide.fa
    
    # Extract the S protein's translated sequence from the GenBank file
    bio view NC_045512:S --trans > protein.fa

    # Extract the S protein's DNA coordinates from the GenBank file
    bio view NC_045512:S --gff > results.gff

## Usage

```{bash, comment=NA}
bio view NC_045512 | head
```

## Command line help
```{bash, comment=NA}
bio view -h
```
