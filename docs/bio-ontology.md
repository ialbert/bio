# `define`: explains terms {#bio-ontology}

The `bio` package provides utility to search gene and sequence ontology.

## Download the databases

Before the first use the databases need to be build or downloaded. If you haven't done so before run this once.

    bio --download

## Database information
 
run
    bio define

prints:

    OntologyDB: total=46,271 gene=43,917 sequence=2,354

There are a total of `47,271` ontology terms out of which `43,917` are gene ontology and `2,529` are sequence ontology terms.

## Define a term 

    bio define exon

prints:

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

## Define term by SO

    bio define SO:0000147

## Define more complex terms

    bio define positive regulation of cell motility

    bio define cellular response to tumor cell

    bio --define intergenic mrna trans splicing

## Building the database

To generate an up to date database use:

    bio define --build

The command above has to be run once (perhaps on a monthly basis) to download the latest data. The efficiency of the process depends on the speed of the hard drive and takes around 30 seconds.

## Showing the term lineage

    bio define exon --lineage

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

Any query that is not matched will be searched for, 
The `-go` flag filters for gene ontology while  `-so` filters for sequence ontology.

Without the `-so` or `-go` flags, it will print out both.

To search for both sequence and gene ontology:

    bio define histone | head 

## Search gene ontology only

    bio define histone --go | head 


To search for sequence ontology:

    bio define histone --so |head




