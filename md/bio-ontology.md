# Taxonomy operations {#bio-taxonomy}

The `bio` package provides utility to visualize NCBI taxonomies.

## Building the taxonomy

Before using the taxonomy related functionality the representation needs to be built:

    bio taxon --download --build

The command above has to be run once (perhaps on a monthly basis) to download and process the latest NCBI taxonomy. The efficiency of the process depends on the speed of the hard drive and takes around 30 minutes.

## Check database


    bio --taxon
    
prints:

    TaxDB: nodes=2,288,072 parents=198,666

There are a total of `2,288,072` nodes (taxonomical entries) out of which `198,666` are nodes that are non-terminal (non-leaf) nodes. For these numbers we see that the vast majority of the taxonomy annotations are for terminal, leaf nodes.

## Searching for taxids

Searches the taxonomy for a word

    bio taxon jawed 

prints:

    # searching taxonomy for: jawed
    clade, Gnathostomata (jawed vertebrates), 7776
    species, Gillichthys mirabilis (long-jawed mudsucker), 8222
    species, Pseudamia amblyuroptera (white-jawed cardinalfish), 1431476
    species, Myctophum brachygnathum (short-jawed lanternfish), 1519985
    species, Oryzias orthognathus (sharp-jawed buntingi), 1645897
    species, Longjawed orbweaver circular virus 1, 2293294
    species, Longjawed orbweaver circular virus 2, 2293295

The search words may use regular expression control characters:

    bio taxon '^jawed'

produces:

    # searching taxonomy for: ^jawed
    clade, Gnathostomata (jawed vertebrates), 7776

## View taxonomy for data 

Once you fetch the data
    
    bio NC_045512 --fetch --rename ncov
        
you can view the descendants:

```{bash, comment=NA}
bio ncov --taxon
```

or view the lineage:

```{bash, comment=NA}
bio ncov --taxon --lineage
```

## View taxonomy by tax id
    
Pass a NCBI taxonomical id to see all the descendants of it:

```{bash, comment=NA}
bio 117565 --taxon | head
```

## View a tax id 

Pass a NCBI taxonomical id to see all the descendants of it:

```{bash, comment=NA}
bio 117565 --taxon | head
```

To print the lineage of a term use:

```{bash, comment=NA}
    bio 564286 --taxon --lineage
```

the lineage may be flattened:

```{bash, comment=NA}
    bio 564286 --taxon --lineage --flat
```
   
## Filter blast results

(TODO) - filters BLAST alignment for species that fall within a taxonomical clade

## List the content of the database:

    bio --taxon --list | head
    
prints:

    no rank, root, 1
    superkingdom, Bacteria (eubacteria), 2
    genus, Azorhizobium, 6
    species, Azorhizobium caulinodans, 7
    species, Buchnera aphidicola, 9
    genus, Cellvibrio, 10
    species, Cellulomonas gilvus, 11
    genus, Dictyoglomus, 13
    species, Dictyoglomus thermophilum, 14
    genus, Methylophilus, 16

Note: this command benefits greatly from using `--preload`.

## Preloading data

For many use cases,  the default behavior is plenty fast and can produce family, genus and species level information in a fraction of a second.

Internally, during operation, the software will query the database for each child node. When selecting a rank where the number of descendant nodes is large (over 10,000 nodes) the run time of the independent queries adds up to a substantial overhead.

For example the command below attempts to render the complete NCBI taxonomic tree with over 2.2 million descendant nodes. When run like so it will take a very long time to produce the output (more than two hours):

    bio 1 --taxon 

The software can operate in a different mode to speed up the process massively by preloading all the data into memory at the cost of imposing a 6 second pre-loading penalty.

    bio 1 --taxon --preload
    
When run with the `--preload` flag the command takes a total of just 11 seconds to generate the same large tree of the entire NCBI taxonomical tree. We don't apply this mode by default because all queries would then take at least 6 seconds, even those that currently finish very quickly.

For queries that take more than 10 seconds to complete (have more than 10,000 descendant nodes) we recommend applying the `--preload` flag.

 