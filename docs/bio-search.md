# bio search: find information {#bio-search}

`bio search` provides command line interface to SRA and MyGene. First install `bio` with:

    pip install bio --upgrade

## Searching SRA (Short Read Archive)

SRA stores information by SRA accession numbers. Once you know the accession use it like so:

    bio search SRR14575325

prints:

    [
        {
            "run_accession": "SRR14575325",
            "sample_accession": "SAMN19241174",
            "first_public": "2021-05-20",
            "country": "",
            "sample_alias": "GSM5320434",
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

## Search BioProjects

    bio search PRJNA257197 > data.json

generates a file with 891 SRR runs.

## How to change JSON into simpler formats

The output of the `bio search` command is in the JSON (JavaScript Object Notation) format.  JSON is not biology specific, instead, you can think of it as a generic, lightweight, human readable approach to represent diverse and hierarchical information.

* https://en.wikipedia.org/wiki/JSON

A tool called `jq` may be used to extract information from JSON data.

    mamba install jq

We use `jq` by forming extraction patterns with the following rules:

* key names are specified with a period `.`
* lists with square brackets `[]`.

We build the extraction command left to right, matching various aspects of it. Since the output is always a list, our first pattern will usually be `.[]` that visits all elements in the list.

For example, the pattern  `.[].run_accession` instructs `jq` to fetch the list for the root key `.[]`, then visit each element and extract the value for the key named `run_accession`. Use the pattern it like so:

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
