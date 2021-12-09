# bio search: find information {#bio-search}

`bio search` provides command line interface to SRA and MyGene. Install `bio` with:

    # Install the package
    pip install bio --upgrade

    # Update the database
    bio --download

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.

## Quick start

The `search` command attempts to provide a unified interface into a multitude of data source:

    # Search nucleotide GenBank
    bio search AF086833

    # Search protein GenBank
    bio search NP_000509

    # Search SRA for a bioproject information
    # Return max of 10 rows in tab delimited format
    bio search  PRJNA257197 --limit 10 -tab

    # Search SRA for a run information
    bio search SRR14575325

    # Search SRR and format output as comma separated values
    bio search SRR14575325 --csv

    # Search MyGene for gene symbol.
    bio search symbol:HBB --limit 1

    # Search MyGene and annotate with Ensembly ids.
    bio search symbol:HBB --species human --fields ensembl

    # Search the NCBI assemblies by a genome build number
    bio search GCA_000002415

    # Search the NCBI assemblies by an organism name, return tab delimited results.
    bio search yoelii -tab

The returned value will be in JSON.  For some searches it may also be comma separated values (`-csv`) or tab separated values (`-tab`)

## How does this tool work?

`bio search` attempts to recognize the search words with bioinformatics signficance (like accession numbers) and attempts to look them up in the most appropriate database. If the search term does not match any known format, it searches the GenBank assemblies for regular expression matches. First install the downloadable database withL

Then start searching for information:

    bio search AF086833

will produce the output:

    [
        {
            "accessionversion": "AF086833.2",
            "assemblyacc": "",
            "assemblygi": "",
            "biomol": "cRNA",
            "biosample": "",
            "caption": "AF086833",
            "completeness": "complete",
            "createdate": "1999/02/10",
            "extra": "gi|10141003|gb|AF086833.2|",
            "flags": "",
            "geneticcode": "1",
            "genome": "",
            "gi": 10141003,
            "moltype": "rna",
            "organism": "Ebola virus - Mayinga, Zaire, 1976",
            "projectid": "0",
            "segsetsize": "",
            "slen": 18959,
            "sourcedb": "insd",
            "strain": "Mayinga",
            "strand": "",
            "subname": "Mayinga|EBOV-May",
            "subtype": "strain|gb_acronym",
            "taxid": 128952,
            "tech": "",
            "term": "10141003",
            "title": "Ebola virus - Mayinga, Zaire, 1976, complete genome",
            "topology": "linear",
            "uid": "10141003",
            "updatedate": "2012/02/13"
        }
    ]

You may also turn the output into comma separated or tab delimited formats:

    bio search AF086833 --csv

produces CSV files:

    AF086833.2,,,cRNA,,AF086833,complete,1999/02/10,gi|10141003|gb|AF086833.2|,,1,,10141003,rna,"Ebola virus - Mayinga, Zaire, 1976",0,,18959,insd,Mayinga,,Mayinga|EBOV-May,strain|gb_acronym,128952,,10141003,"Ebola virus - Mayinga, Zaire, 1976, complete genome",linear,10141003,2012/02/13

whereas to produce tab delimited formats use:

    bio search AF086833 --tab

prints:

    AF086833.2  cRNA  AF086833  complete  1999/02/10  gi|10141003|gb|AF086833.2|  1  10141003  rna  Ebola  virus  -  Mayinga,  Zaire,  1976  0  18959  insd  Mayinga  Mayinga|EBOV-May  strain|gb_acronym  128952  10141003  Ebola  virus  -  Mayinga,  Zaire,  1976,  complete  genome  linear  10141003  2012/02/13


## Search of SRR numbers

    bio search SRR14575325

produces:

    [
        {
            "run_accession": "SRR14575325",
            "sample_accession": "SAMN19241174",
            "first_public": "2021-05-20",
            "country": "",
            "sample_alias": "GSM5320434",
            "fastq_bytes": "612817621",
            "read_count": "19439295",
            "library_name": "",
            "library_strategy": "miRNA-Seq",
            "library_source": "TRANSCRIPTOMIC",
            "library_layout": "SINGLE",
            "instrument_platform": "ILLUMINA",
            "instrument_model": "Illumina HiSeq 2000",
            "study_title": "Identification of 5'isomiR in HCC patients.",
            "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/025/SRR14575325/SRR14575325.fastq.gz"
        }
    ]

To search the MyGene interface add a type to the word:

    bio search symbol:HBB --limit 1

it will print:

    [
        {
            "name": "hemoglobin subunit beta",
            "refseq": {
                "genomic": [
                    "NC_000011.10",
                    "NG_000007.3",
                    "NG_059281.1"
                ],
                "protein": "NP_000509.1",
                "rna": "NM_000518.5",
                "translation": {
                    "protein": "NP_000509.1",
                    "rna": "NM_000518.5"
                }
            },
            "symbol": "HBB",
            "taxid": 9606,
            "taxname": "Homo sapiens"
        }
    ]
    #  showing 1 out of 30 results.

If the query word does not seem to match any known format a regular expression search will applied on the GenBank assembly reports:

    bio search plasmodium

prints:

    [
        {
            "assembly_accession": "GCA_000002415.2",
            "bioproject": "PRJNA150",
            "biosample": "SAMN02953638",
            "wgs_master": "AAKM00000000.1",
            "refseq_category": "representative genome",
            "taxid": "5855",
            "species_taxid": "5855",
            "organism_name": "Plasmodium vivax",
            "infraspecific_name": "",
            "isolate": "Salvador I",
            "version_status": "latest",
            "assembly_level": "Chromosome",
            "release_type": "Minor",
            "genome_rep": "Full",
            "seq_rel_date": "2009/05/06",
            "asm_name": "ASM241v2",
            "submitter": "TIGR",
            "gbrs_paired_asm": "GCF_000002415.2",
            "paired_asm_comp": "different",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/415/GCA_000002415.2_ASM241v2",
            "excluded_from_refseq": "",
            "relation_to_type_materialasm_not_live_date": ""
        },
        ...
    ]

## Other outputs

The output may be set to comma separated values or tab delimited values

    bio search plasmodium --csv

    bio search plasmodium --tab

## How to transform JSON into other formats

Our tabular formats may be insufficient when to content is nested. In those cases you will need to transform the JSON into other formats.

The output of the `bio search` command is produced the so called JSON (JavaScript Object Notation) format.  JSON is not biology specific, instead, you can think of it as a generic, lightweight, human readable approach to represent diverse and hierarchical information.

* https://en.wikipedia.org/wiki/JSON

A tool called `jq` may be used to extract information from JSON data in various ways.

    mamba install jq

We use `jq` by forming extraction patterns with the following rules:

* key names are specified with a period `.`
* lists with square brackets `[]`.

We build the extraction command left to right, matching various aspects of it. Since the output is always a list, our first pattern will usually be `.[]` that visits all elements in the list.

For example, the pattern  `.[].run_accession` instructs `jq` to fetch the list for the root key `.[]`, then visit each element and extract the value for the key named `run_accession`. Use the pattern it like so:

    # Produce a JSON output.
    bio search PRJNA257197 > data.json

    # Extract information from the JSON.
    cat data.json | jq -r .[].run_accession | head -5

prints the accession numbers:

    SRR1553416
    SRR1553417
    SRR1553418
    SRR1553419
    SRR1553420

To extract two colums `run_accession` and `read_count` into a table:

    cat data.json | jq -r '.[]|[.run_accession,.read_count]|@tsv' | head -5

prints:

    SRR1553416	1014800
    SRR1553417	1865865
    SRR1553418	18323
    SRR1553419	281432
    SRR1553420	223206

When I don't use `jq` for a while, I forget how it works. I found that the example above helps remind me and is a good starting point to build other patterns.

## How to learn `jq`

The way `jq` works may not be everyone's cup of tea. On the other hand you can do wondrous things with `jq` if you put in some effort to learn it. Stick with it! Once you get it, making patterns becomes that puzzle-solving kind of fun.

There are numerous resources on learning `jq`. Search and find one that suits you. I found the page below quite useful:

* https://earthly.dev/blog/jq-select/

## The MyGene interface

If the search terms do not match the SRA number format the search proceeds with querying the MyGene interface. The `mygene` service was published in [High-performance web services for querying gene and variant annotation, Genome Research, 2016][mygene] and is an amazingly potent query interface that focuses on returning data rather than web pages.

Our search is a command line convenience function that stands on the shoulders of the giants that created the [mygene service][mygene] in the first place.

[mygene]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9

## Usage example

    bio search HBB --limit 1

will print:

    [
        {
            "name": "hemoglobin subunit beta",
            "symbol": "HBB",
            "taxid": 9606,
            "taxname": "Homo sapiens"
        }
    ]
    #  showing 1 out of 148 results.

The same search for `HBB` will return a lot more data as we add fields to it. Lets focus the search on symbol only `symbol:HBB`  human as species, and request Ensembl annotations:

    bio search symbol:HBB --species human --fields ensembl

the query will now produce a much larger dataset:

    [
        {
            "ensembl": {
                "gene": "ENSG00000244734",
                "protein": [
                    "ENSP00000333994",
                    "ENSP00000369671",
                    "ENSP00000488004",
                    "ENSP00000494175",
                    "ENSP00000496200"
                ],
                "transcript": [
                    "ENST00000335295",
                    "ENST00000380315",
                    "ENST00000475226",
                    "ENST00000485743",
                    "ENST00000633227",
                    "ENST00000647020"
                ],
                "translation": [
                    {
                        "protein": "ENSP00000333994",
                        "rna": "ENST00000335295"
                    },
                    {
                        "protein": "ENSP00000494175",
                        "rna": "ENST00000647020"
                    },
                    {
                        "protein": "ENSP00000488004",
                        "rna": "ENST00000633227"
                    },
                    {
                        "protein": "ENSP00000496200",
                        "rna": "ENST00000485743"
                    },
                    {
                        "protein": "ENSP00000369671",
                        "rna": "ENST00000380315"
                    }
                ],
                "type_of_gene": "protein_coding"
            },
            "name": "hemoglobin subunit beta",
            "symbol": "HBB",
            "taxid": 9606,
            "taxname": "Homo sapiens"
        }
    ]

Below is a more complicated construct with `jq` that pairs up fields into a tab-delimited table.

    cat data.json | jq -r '.[].ensembl.translation[]|[.protein,.rna]| @tsv'

    ENSP00000333994    ENST00000335295
    ENSP00000494175    ENST00000647020
    ENSP00000488004    ENST00000633227
    ENSP00000496200    ENST00000485743
    ENSP00000369671    ENST00000380315


## Additional examples

Often examples are more useful than words

    # Limiting searches
    bio search HBB --limit 10

    # Adding summary
    bio search HBB -f summary --limit 1 --species human

    # Adding additional fields.
    bio search symbol:HBB --species human --fields Ensembl,RefSeq

    # Gene products with the function in human
    bio search go:0000307 --species human

    # Text search for a concept.
    bio search insulin

    # Genes with insulin in the summary then also print the summary field
    bio search summary:insulin -f summary

    # Find aliases, produces the same output
    bio search symbol:XRCC2 -f alias --species human
    bio search alias:SPGF50 -f alias --species human

## Field documention

See the fields listed at API documentation site:

* https://docs.mygene.info/en/latest/index.html

* https://docs.mygene.info/en/latest/doc/query_service.html#available_fields
