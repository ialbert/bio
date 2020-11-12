# Taxonomy {#taxon}

The `bio` provides utility to visualize NCBI taxonomies.

## Building the taxonomy

Before using the taxonomy related functionality the representation needs to be built:

    bio taxon --download --build

The command above will download and process the NCBI taxonomy to prepare it for fast access.

## Using taxonomies

Pass an NCBI taxid to view the child nodes:

```{bash, comment=NA}
bio taxon 117565 | head
```

## Preloading into memory

When the output of the taxonomy is very large, for example the command below
will render the complete NCBI tree with over 2 million entries

    bio taxon 1

the output is starts streaming immediately but the database is hit for each node, hence slowing
down the data generation.

It is possible to preload all the data in memory to speed up the process at the cost of imposing a 6 second loading penalty.

    bio taxon 1 --preload
    
The command above takes a total of just 10 seconds to generate a large tree with the entire 
NCBI taxonomical tree.    


## Searching for taxids

(TODO) - produces a taxid when searching for a word

## Render lineages

(TODO) - renders the complete lineage of a term


## Filter blast results

(TODO) - filters BLAST alignment for species that fall within a taxonomical clade
    