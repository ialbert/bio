# bio mygene: get gene information {#bio-mygene}

The `bio mygene` command provides command-line based access to the https://mygene.info/ query interface.

The `mygene` service was published in [High-performance web services for querying gene and variant annotation, Genome Research, 2016][mygene] and is an amazingly potent query interface that focuses on returning data rather than web pages.

Note that `bio mygene` is merely a command line convenience function that stands on the shoulders of the giants that created the [mygene service][mygene] in the first place.

[mygene]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9

## Rationale

We have implemented the `bio mygene` command to facilitate quick explorations of gene information. For more information on data sources and representations, consult [The Biostar Handbook][book] chapters on [Biological Data Sources][datasource]. To install `bio` use:

[datasource]: https://www.biostarhandbook.com/biological-data-sources.html
[book]: https://www.biostarhandbook.com

    pip install bio --upgrade

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.

## Usage example

    bio mygene HBB --limit 1

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

The same search for `HBB` will return a lot more data as we add fields to it. Let's focus the search on symbol only `symbol:HBB` human as species, and request Ensembl annotations:

    bio mygene symbol:HBB --species human --fields ensembl

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

The output of the `bio mygene` command is in the JSON (JavaScript Object Notation) format. JSON is not biology-specific, instead, you can think of it as a generic, lightweight, human-readable approach to represent diverse and hierarchical information.

- https://en.wikipedia.org/wiki/JSON

## How to change JSON into simpler formats

The tool called `jq` may be used to extract information from JSON data.

    mamba install jq

to use `jq` first store the output in a file:

    bio mygene symbol:HBB --species human --fields ensembl > data.json

Then practice the building extraction patterns with the following rules:

- key names are specified with a period `.`
- lists with square brackets `[]`.

We build the extraction command left to right, matching various aspects of it. Since the output is always a list, our first pattern will usually be `.[]` that visits all elements in the list.

For example, the pattern `.[].symbol` instructs `jq` to fetch the list for the root key `.[]`, then visit each element and extract the value for the key named `symbol`. Use the pattern like so:

    cat data.json | jq -r '.[].symbol'

prints:

    HBB

Here is a more complicated pattern: `.[].ensembl.protein[]` it will visit each element of the root list, extract the value for the `ensemble` key, then the value for the `protein` key and then prints all the resulting elements of that list:

    cat data.json | jq -r '.[].ensembl.protein[]'

    ENSP00000333994
    ENSP00000369671
    ENSP00000488004
    ENSP00000494175
    ENSP00000496200

## Advanced usage

`jq` is quite powerful. Below is a more complicated construct that pairs up fields into a tab-delimited table.

    cat data.json | jq -r '.[].ensembl.translation[]|[.protein,.rna]| @tsv'

    ENSP00000333994    ENST00000335295
    ENSP00000494175    ENST00000647020
    ENSP00000488004    ENST00000633227
    ENSP00000496200    ENST00000485743
    ENSP00000369671    ENST00000380315

Well, yeah ... not everyone's cup of tea, suffice to say you can do wondrous things with `jq` if you put in some effort to learn it.

Stick with it! Once you get it, making patterns becomes that puzzle-solving kind of fun.

## How to learn `jq`

There are numerous resources on learning `jq`. Search and find one that suits you. I found the page below quite useful:

- https://earthly.dev/blog/jq-select/

## Additional examples

Often examples are more useful than words

    # Limiting searches
    bio mygene HBB --limit 10

    # Adding summary
    bio mygene HBB -f summary --limit 1 --species human

    # Adding additional fields.
    bio mygene symbol:HBB --species human --fields Ensembl,RefSeq

    # Gene products with function in human
    bio mygene go:0000307 --species human

    # Text search for a concept.
    bio mygene insulin

    # Genes with insulin in the summary then also print the summary field
    bio mygene summary:insulin -f summary

    # Find aliases, produces the same output
    bio mygene symbol:XRCC2 -f alias --species human
    bio mygene alias:SPGF50 -f alias --species human

## Field documentation

See the fields listed at the API documentation sites:

- https://docs.mygene.info/en/latest/index.html

- https://docs.mygene.info/en/latest/doc/query_service.html#available_fields
