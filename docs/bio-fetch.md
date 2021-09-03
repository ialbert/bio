# `bio fetch`: get data {#bio-fetch}

`bio` can fetch data from GenBank or from Ensembl.

## Fetch data from GenBank

You may use multiple accession numbers:

    bio fetch NC_045512 MN996532 > genomes.gb

You may pipe accession numbers into the tool

    echo NC_045512 | bio fetch > genomes.gb

the `fetch` task is not a replacement for other means of accessing NCBI, notably `entrez-direct`. Instead, think of it as a convenience function that simplifies a few common use cases.

For more advance command line data access options see Entrez Direct:

* https://www.ncbi.nlm.nih.gov/books/NBK179288/

## Fetch data from Ensembl

`bio fetch` recognzies Ensemble gene and transcript names (ENSG, ENST) and will automatically connect to the Ensembl REST API:

    bio fetch ENSG00000157764 | head

    # Transcript data un genomic context
    bio fetch ENST00000288602  | head

    # Transcript data as CDNA
    bio fetch ENST00000288602 --type cdna | head

    # Transcript data as CDS
    bio fetch ENST00000288602 --type cds | head

    # Transcript data as protein
    bio fetch ENST00000288602 --type cds | head


For more information see the Ensembl REST API:

* https://rest.ensembl.org/
