# bio fetch: download data {#bio-fetch}

We have implemented the `bio fetch` command to facilitate data download from GenBank, SRA and Ensembl. Some of the fetch commands  build upon or directly rely on functionality already present in other tools such as Entrez Direct:

* [Entrez Direct: E-utilities on the Unix Command Line][entrez-direct]

If you need additional functionality beyond of what `bio fetch` offers see the end of the page for a list of alternatives.

[entrez-direct]: https://www.ncbi.nlm.nih.gov/books/NBK179288/

For more information on data sources and representations, consult [The Biostar Handbook][book] chapters on [Biological Data Sources][datasource]. To install `bio` use:

[datasource]: https://www.biostarhandbook.com/biological-data-sources.html
[book]: https://www.biostarhandbook.com

    pip install bio --upgrade
    bio --download

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.

## Rationale

`bio fetch` was written to provide a simpler data access. Since in most cases the accession numbers already uniquely define the resource we can implement a tool that bypasses the additional configuration that official data sources demand. Note how much simpler the following is:

    # Obtains nucleotidedata from GenBank
    bio fetch NC_045512

    # Obtains protein data from GenBank
    bio fetch YP_009724390

    # Obtains FASTA sequence files from Ensembl
    bio fetch ENSG00000157764

    # Obtain the transcript data as CDS
    bio fetch ENST00000288602 --type cds | head

    # Obtains SRR run information from the Short Read Archive
    bio fetch SRR1972976

## Fetch data from GenBank

Get Genbank nucleotides by accession number

	bio fetch NC_045512 | head

`fetch` automatically recognizes protein ids and connects to protein database, no further configuration is needed:

	bio fetch YP_009724390 | head

You may also list multiple accession numbers:

    bio fetch NC_045512 MN996532 > genomes.gb

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

And the [enaBrowserTools: interface with the ENA web services to download data from ENA][ena]

## Fetch run information from BioProject id

[PRJNA257197]:https://www.ncbi.nlm.nih.gov/bioproject/PRJNA257197/

The following command produces a comma separated run information associated with a bioproject id, for example take [BioProject:PRJNA257197]:

    bio fetch PRJNA257197 --limit 1

will produce the comma separated output:

    Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
    SRR1972976,2015-04-14 13:53:37,2015-04-14 13:48:38,8345287,1685747974,8345287,202,997,,https://sra-downloadb.st-va.ncbi.nlm.nih.gov/sos2/sra-pub-run-6/SRR1972976/SRR1972976.1,SRX994253,W220.0.l1,RNA-Seq,cDNA,TRANSCRIPTOMIC,PAIRED,0,0,ILLUMINA,Illumina HiSeq 2500,SRP045416,PRJNA257197,2,257197,SRS908478,SAMN03253746,simple,186538,Zaire ebolavirus,W220.0,,,,,,,no,,,,,BI,SRA178666,,public,6A26DBAB1096535FCB94FCE9E1AE8AD8,FB20A0391119E532EA03F374A16EB508

Use `csvcut` to select columns of interest:


    bio fetch PRJNA257197 --limit 1 | csvcut -c Run,Platform,spots,avgLength

prints:

    Run,Platform,spots,avgLength
    SRR1972976,ILLUMINA,8345287,202

The `bio fetch PRJNA257197 --limit 1` command is equivalent to running `entrez-direct` construct:

    esearch -db sra -query PRJNA257197 | efetch -stop 1 -format runinfo

## Fetch SRA run information

[sra]:https://www.ncbi.nlm.nih.gov/sra

The following command retrievs JSON data that describes sequencing data deposited at the [Short Read Archive][sra]

    bio fetch SRR1972976 | csvcut -c Run,ScientificName,TaxID

produces:

    Run,ScientificName,TaxID
    SRR1972976,Zaire ebolavirus,186538

the command `bio fetch SRR1972976` is equivalent to running

    efetch -db sra -id SRR1972976 -format runinfo

## Tools with similar utility {#fetch_similar}

`bio fetch` is primarily a convenience function that simplifies the use of entrez in certain simple and well defined cases. For all other cases please consult the documentation for Entrez Direct:

* [Entrez Direct: E-utilities on the Unix Command Line][entrez-direct]


See also the related tools that may have expanded functionality:

* [ffq: Fetch run information from the European Nucleotide Archive (ENA)][ffq]
* [enaBrowserTools: interface with the ENA web services to download data from ENA][ena]

[ffq]: https://github.com/pachterlab/ffq
[ena]: https://github.com/enasequence/enaBrowserTools
