# bio explain: show definitions {#bio-explain}

We implemented the `bio explain` command to facilitate the quick explorations of the Gene and Sequence Ontologies.

For more information on ontologies, consult [The Biostar Handbook][book] chapter [Sequence and Gene Ontology: What do the words mean?][ontology]. To install `bio` use:

[ontology]: https://www.biostarhandbook.com/what-do-the-words-mean.html
[book]: https://www.biostarhandbook.com

    pip install bio --upgrade
    bio --download

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.

## Database information

    bio explain

Prints the contents of the database:

    # Content: 43,878 gene ontology terms; 2,354 sequence ontology terms

As we can see the database currently contains `43,878` gene ontology and `2,354 `  sequence ontology terms.

## Explain a term

    bio explain exon

prints that the term is part of the Sequence Ontology (SO) with number `SO:0000147` and is defined as:

    ## exon (SO:0000147)

    A region of the transcript sequence within a gene which is not removed from the
    primary RNA transcript by RNA splicing.

    Parents:
    - transcript_region

    Children:
    - coding_exon
    - noncoding_exon
    - interior_exon
    - decayed_exon (non_functional_homolog_of)
    - pseudogenic_exon (non_functional_homolog_of)
    - exon_region (part_of)
    - exon_of_single_exon_gene

You can now query the parent or child relationships:

    bio explain transcript_region

or

    bio explain coding_exon

## Explain term by identifier

If you know the identifier then you may use it directly `bio explain SO:0000147` or `bio explain GO:2000147`

## Explain more complex terms

All words have to match exactly; thus, try simpler terms before going full length.

    bio explain regulation

    bio explain positive regulation

    bio explain positive regulation of x

    bio explain positive regulation of xanthophore differentiation

What is xantophore?

    bio explain xanthophore

prints:

    ## xanthophore (GO:0031633)

    A chromatophore containing yellow pigment.

    Parents:
    - plasma membrane-derived chromatophore

## Showing the term lineage

    bio explain exon --lineage

prints:

    SO:0000110  sequence_feature
      SO:0000001  region
        SO:0001411  biological_region
          SO:0000833  transcript_region

            ## exon (SO:0000147)

            A region of the transcript sequence within a gene which is not removed from the
            primary RNA transcript by RNA splicing.


            Children:
            - coding_exon
            - noncoding_exon
            - interior_exon
            - decayed_exon (non_functional_homolog_of)
            - pseudogenic_exon (non_functional_homolog_of)
            - exon_region (part_of)
            - exon_of_single_exon_gene


## Searching the database

Any query that is not matched will be searched for. The `-go` flag restricts tje search for gene ontology while  `-so` restricts the search for sequence ontology. To search for both sequence and gene ontology:

    # Searces both ontologies
    bio explain histone | head

    # Search gene ontology only
    bio explain histone --go | head

    # Search gene ontology only
    bio explain histone --so |head

## Build the newest version locally

You may build a database with the newest data:

    bio explain --build

The command will download the most up-to-date ontology data and build a new database. The process might take about half an hour - through the completion speed depends on the hard drive write speed, solid-state drives finish much faster.
