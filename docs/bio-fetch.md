# `bio fetch`: get GenBank data {#bio-fetch}

`bio` can fetch data in GeneBank format based on accession numbers

You may separate accession numbers with commas

    bio fetch NC_045512,MN996532 > genomes.gb

You may separate accession numbers with spaces

    bio fetch NC_045512 MN996532 > genomes.gb

You may pipe accession numbers into the tool

    echo NC_045512 | bio fetch > genomes.gb

the `fetch` task is not a replacement for other means of accessing NCBI, notably `entrez-direct`.

Instead think of it as a convenience function that simplifies the most common usecase.


## Fetch data from Ensembl

`bio fetch` recognzies Ensemble gene and transcript names (ENSG, ENST) and will automatically connect to the Ensembl REST API:

    bio fetch ENSG00000157764 | head

    bio fetch ENST00000288602.11 | head


Ensembl REST API:

* https://rest.ensembl.org/
