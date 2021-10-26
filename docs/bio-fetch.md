# bio fetch: download data {#bio-fetch}

We have implemented the `bio fetch` command to facilitate data download from GenBank, Ensembl and other sources.  If you need additional functionality beyond of what `bio fetch` offers see the end of the page for a list of alternatives. Install `bio` with:

    pip install bio --upgrade

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.

## Rationale

`bio fetch` was written to provide simpler data access. As in most cases the accession numbers uniquely define the data type `bio fetch` can automatically resolve the correct destination and will download data in the most appropriate format. Examples:

    # Obtains nucleotidedata from GenBank
    bio fetch NC_045512

    # Obtains protein data from GenBank
    bio fetch YP_009724390

    # Obtains SRR run information from the Short Read Archive
    bio fetch SRR1972976

    # Obtains FASTA sequence files from Ensembl
    bio fetch ENSG00000157764

    # Obtain the coding sequence segments of a transcript
    bio fetch ENST00000288602 --type cds | head

Read on for more details.

## Fetch data from GenBank

To get a Genbank nucleotides by accession number run:

    # Obtains nucleotidedata from GenBank
	bio fetch NC_045512 | head

The default format is GenBank. You may also list multiple accession numbers at once:

    # Obtains mulitple entries from GenBank
    bio fetch NC_045512 MN996532 > genomes.gb

you may also fetch FASTA and GFF data directly from GenBank. The command

    bio fetch NC_045512 --format fasta  | head -3

returns a FASTA file:

    >NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA
    CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC

whereas:

    bio fetch NC_045512 --format gff  | head -3

returns a GFF file:

    ##sequence-region NC_045512.2 1 29903
    ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049
    NC_045512.2	RefSeq	region	1	29903	.	+	.	ID=NC_045512.2:1..29903;Dbxref=taxon:2697049;collection-date=Dec-2019;country=China;gb-acronym=SARS-CoV-2;gbkey=Src;genome=genomic;isolate=Wuhan-Hu-1;mol_type=genomic RNA;nat-host=Homo sapiens;old-name=Wuhan seafood market pneumonia virus

## Fetching genome assemblies

Genome assemblies at NCBI are also accessible with:

    # Produces a GenBank file
    bio search GCA_000006155

    # Produces a FASTA file
    bio search GCA_000006155 --format fasta

    # Produces a GFF file.
    bio search GCA_000006155 --format gff

## Fetch data from Ensembl

`bio fetch` recognizes Ensemble accessions (ENSG, ENST) and will automatically connect to the Ensembl REST API:

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

These commands fetch data in FASTA formats. For more information see the Ensembl REST API:

* https://rest.ensembl.org/

And the [enaBrowserTools: interface with the ENA web services to download data from ENA][ena]


## Tools with similar utility {#fetch_similar}

`bio fetch` is primarily a convenience function that simplifies the use of Entrez and ENA in certain simple and well defined cases.

See also the related tools that may have expanded functionality:

* [Entrez Direct: E-utilities on the Unix Command Line][entrez-direct]
* [enaBrowserTools: interface with the ENA web services to download data from ENA][ena]
* [ffq: Fetch run information from the European Nucleotide Archive (ENA)][ffq]


From https://www.biostars.org/p/9492531/#9492581 use the [ENA portal API](https://www.ebi.ac.uk/ena/portal/api/)

    # Get information
    curl -X GET "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA374918&fields=run_accession,library_name,library_layout,fastq_ftp&result=read_run"

    # What fields are acceptable
    curl -X GET "https://www.ebi.ac.uk/ena/portal/api/returnFields?dataPortal=ena&result=read_run"

[ffq]: https://github.com/pachterlab/ffq
[ena]: https://github.com/enasequence/enaBrowserTools
[entrez-direct]: https://www.ncbi.nlm.nih.gov/books/NBK179288/
