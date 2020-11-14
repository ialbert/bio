# NCBI Taxonomy {#taxon}

The `bio` package provides utility to visualize NCBI taxonomies.

## Building the taxonomy

Before using the taxonomy related functionality the representation needs to be built:

    bio taxon --download --build

The command above will download and process the NCBI taxonomy to prepare it for fast access.

## Quick database check

    bio taxon
    
prints:

    TaxDB: nodes=2,288,072 parents=198,666

There are a total of `2,288,072` nodes (taxonomical entries) out of which `198,666` are nodes that are non-terminal (leaf) nodes that have
children nodes attached to them. Above we see that the vast majority of the taxonomy annotates leaf nodes.

## View a tax id

Pass a NCBI taxonomical id to see all the descendants of it:

```{bash, comment=NA}
bio taxon 117565 | head
```

## Print lineages

To print the lineage of a term use:

```{bash, comment=NA}
    bio taxon 564286 --lineage
```

the lineage may be flattened:

```{bash, comment=NA}
    bio taxon 564286 --lineage --flat
```

it is also possible to pipe data into the tool:

    cat taxids | bio taxon --lineage 
    
## Searching for taxids

Searches the taxonomy for a word

    bio taxon jaw | head

prints:

    # searching taxonomy for: jaw
    clade, Gnathostomata (jawed vertebrates), 7776
    species, Gillichthys mirabilis (long-jawed mudsucker), 8222
    family, Oplegnathidae (knifejaws), 30858
    species, Gonostoma atlanticum (Atlantic fangjaw), 48456
    species, Sigmops gracilis (slender fangjaw), 48457
    species, Coregonus zenithicus (shortjaw cisco), 56554
    species, Malacosteus niger (stoplight loosejaw), 76143
    species, Gillichthys seta (shortjaw mudsucker), 79683
    species, Galaxias postvectis (shortjaw kokopu), 89561

The search words may use regular expression control characters:

    bio taxon '^jawed'

produces:

    # searching taxonomy for: ^jawed
    clade, Gnathostomata (jawed vertebrates), 7776

## Filter blast results

(TODO) - filters BLAST alignment for species that fall within a taxonomical clade

## List the content of the database:

    bio taxon --list | head
    
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

For many usecases,  the default behavior is plenty fast and can produce family, genus and species level information in a fraction of a second.

Internally, during operation the software will query the database for each child node. When selecting a rank where the number of descendant nodes happens to be large (over 10,000) the independent queries will  add up to a substantial overhead.

For example the command below attempts to render the complete NCBI taxonomic tree that has over 2.2 million descendant nodes. When run like so it will take a long time to produce the output (more than two hours):

    bio taxon 1 

The software can operate in a different mode to speed up the process massively by preloading all the data into memory at the cost of imposing a 6 second loading penalty.

    bio taxon 1 --preload
    
When run with the `--preload` flag the command takes a total of just 11 seconds to generate the same large tree of the entire NCBI taxonomical tree. We don't apply this mode by default because all queries would then take at least 6 seconds, even those that currently finish very quickly.

For queries that take more than 10 seconds to complete (have more than 10,000 descendant nodes) we recommend applying the `--preload` flag.

 