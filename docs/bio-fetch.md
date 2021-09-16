# bio fetch: download data {#bio-fetch}

We have implemented the `bio fetch` command to facilitate data download from GenBank and Ensembl.

For more information on data sources and representations, consult [The Biostar Handbook][book] chapters on [Biological Data Sources][datasource]. To install `bio` use:

[datasource]: https://www.biostarhandbook.com/biological-data-sources.html
[book]: https://www.biostarhandbook.com

    pip install bio --upgrade
    bio --download

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.


## Fetch data from GenBank

Get Genbank nucleotides by accession number

	bio fetch NC_045512 | head

Automatically recognizes protein ids and connects to protein database, no further parameters needed:

	bio fetch YP_009724390 | head

You may also list multiple accession numbers:

    bio fetch NC_045512 MN996532 > genomes.gb

The input accession numbers may be stored in a file. If `acc.txt` contains:

	NC_045512
	MN996532

You may pipe the above accession numbers as standard input:

    catt acc.txt | bio fetch > genomes.gb

the `fetch` task is not a replacement for other means of accessing NCBI, notably `entrez-direct`. Instead, think of it as a convenience function that simplifies a few common use cases.

For more advanced command line data access options to NCBI see Entrez Direct

[bh-entrez]: https://www.biostarhandbook.com/automating-access-to-ncbi.html
[entrez-direct]: https://www.ncbi.nlm.nih.gov/books/NBK179288/

* [Biostar Handbook: Automating access to NCBI][bh-entrez]
* [Entrez Direct: E-utilities on the Unix Command Line][entrez-direct]

## Fetch data from Ensembl

`bio fetch` recognizes Ensemble gene and transcript names (ENSG, ENST) and will automatically connect to the Ensembl REST API:

	# Ensenble gene
    bio fetch ENSG00000157764 | head

    # Transcript data in genomic context
    bio fetch ENST00000288602  | head

    # Transcript data as CDNA
    bio fetch ENST00000288602 --type cdna | head

    # Transcript data as CDS
    bio fetch ENST00000288602 --type cds | head

    # Transcript data as protein
    bio fetch ENST00000288602 --type protein | head

For more information see the Ensembl REST API:

* https://rest.ensembl.org/
