# Ontology operations {#bio-ontology}

The `bio` package provides utility to visualize gene and sequence ontology's.

## Building the database
Before using the ontology related functionality the representation needs to be built:

    bio --define --build

The command above has to be run once (perhaps on a monthly basis) to download the latest data. The efficiency of the process depends on the speed of the hard drive and takes around one minute.

## Check database


    bio --define
    
prints:

    OntologyDB: total=49,747, gene=47,218 sequence=2,529

There are a total of `47,218` ontology terms out of which `47,218` are gene and `2,529` are sequence.

## Define a term 

    bio --define exon
    
prints:

    SO:0000147	exon	"A region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing." [SO:ke]
    SO:0000198	noncoding_exon
    SO:0000201	interior_exon
    SO:0005845	exon_of_single_exon_gene
    SO:0000195	coding_exon
    
    # Define by the SO or GO ID's
    bio --define SO:0000147
    

    SO:0000147	exon	"A region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing." [SO:ke]
    SO:0000198	noncoding_exon
    SO:0000201	interior_exon
    SO:0005845	exon_of_single_exon_gene
    SO:0000195	coding_exon

The first line is the ontological term that matches, with each subsequent line being a child of the first line.

## Showing the term lineage

    bio --define exon --lineage

prints:

    SO:0000110 sequence_feature
        SO:0000001 region
            SO:0001411 biological_region
                SO:0000833 transcript_region
                    SO:0000147 exon
    
    
    # Show using the GO or SO ID.
    bio --define SO:0000147 --lineage
    
prints:

    SO:0000110 sequence_feature
        SO:0000001 region
            SO:0001411 biological_region
                SO:0000833 transcript_region
                    SO:0000147 exon
    

## Searching the database

Any query that is not matched will be searched for, 
The `-go` flag filters for gene ontology while  `-so` filters for sequence ontology.

Without the `-so` or `-go` flags, it will print out both matches.

To search for both sequence and gene ontology:
    
    # Search genes
    bio --define histone 
    
    
prints:

    ...
    ...
    ...
    GO:2001254 negative regulation of histone h3-k36 trimethylation
    GO:2001255 positive regulation of histone h3-k36 trimethylation
    SO:0002143 histone_2b_acetylation_site
    SO:0002144 histone_2az_acetylation_site
    
To search for gene ontology:
    
    # Search genes
    bio --define histone -go 
    
    
prints:

    ...
    ...
    ...
    GO:2001254 negative regulation of histone h3-k36 trimethylation
    GO:2001255 positive regulation of histone h3-k36 trimethylation


To search for sequence ontology:

    # Search sequences
    bio --define histone -so 

prints :

    ...
    ...
    ...
    SO:0002143 histone_2b_acetylation_site
    SO:0002144 histone_2az_acetylation_site
    
This searches the gene, if given `-go`, or sequence, if given `-so`, ontology for `histone` then prints matching results.
Without `-so` or `-go`, it will print out both matches.

