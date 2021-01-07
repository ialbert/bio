# Ontology operations {#bio-ontology}

The `bio` package provides utility to search gene and sequence ontology.

## Building the database
Before using the ontology related functionality the representation needs to be built:

    bio --define --build

The command above has to be run once (perhaps on a monthly basis) to download the latest data. The efficiency of the process depends on the speed of the hard drive and takes around 30 seconds.

## Check database
```{bash, comment=NA}
# Check the database
bio --define
```

There are a total of `47,218` ontology terms out of which `47,218` are gene and `2,529` are sequence.

## Define a term 
```{bash, comment=NA}
# Define the term
bio --define exon
```
```{bash, comment=NA}
# Define term by SO id
bio --define SO:0000147
```
    
The first line is the ontological term that matches, with each subsequent line being a child of the first one.


    bio --define positive regulation of cell motility

    bio --define cellular response to tumor cell

    bio --define intergenic mrna trans splicing

## Showing the term lineage

```{bash, comment=NA}
# show term lineage 
bio --define exon --lineage
```

```{bash, comment=NA}
# Show term lineage by SO name
bio --define SO:0000147 --lineage
```


## Searching the database

Any query that is not matched will be searched for, 
The `-go` flag filters for gene ontology while  `-so` filters for sequence ontology.

Without the `-so` or `-go` flags, it will print out both.

To search for both sequence and gene ontology:

```{bash, comment=NA}
# Search by a keyboard
bio histone --define | head 
```


    
To search for gene ontology:
    
```{bash, comment=NA}
# Search by a keyboard
bio histone --define -go |head 
```


To search for sequence ontology:

```{bash, comment=NA}
# Search by a keyboard
bio histone --define -so |head
```

## Preloading data

For many use cases,  the default behavior is plenty fast and can produce family, genus and species level information in a fraction of a second.

Internally, during operation, the software will query the database for each child node. When selecting a rank where the number of descendant nodes is large (over 10,000 nodes) the run time of the independent queries adds up to a substantial overhead.


When run like so it will around 6 seconds:

    bio regulation --define

The software can operate in a different mode to speed up the process massively by preloading all the data into memory at the cost of imposing a 1 second pre-loading penalty.

    bio regulation --define --preload
    
When run with the `--preload` flag the command takes less than a 2 seconds to generate the same result. 
We don't apply this mode by default because all queries would then take at least 1 second, even those that currently finish very quickly.

For queries that take more than 1 second to complete we recommend applying the `--preload` flag.





