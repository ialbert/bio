# NCBI Taxonomy {#taxon}

The `bio` package provides utility to visualize NCBI taxonomies.

## Building the taxonomy

Before using the taxonomy related functionality the representation needs to be built:

    bio taxon --download --build

The command above will download and process the NCBI taxonomy to prepare it for fast access.

## View a tax id

Pass an NCBI taxonomical id to see all the descendants of it:

```{bash, comment=NA}
bio taxon 117565 | head
```

prints:

```
class, Myxini, 117565
   order, Myxiniformes, 7761
      family, Myxinidae (hagfishes), 7762
         subfamily, Eptatretinae, 30309
            genus, Eptatretus, 7763
               species, Eptatretus burgeri (inshore hagfish), 7764
               species, Eptatretus stoutii (Pacific hagfish), 7765
               species, Eptatretus okinoseanus, 7767
               species, Eptatretus atami, 50612
               species, Eptatretus cirrhatus (broadgilled hagfish), 78394
```

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

The search words may use regular expression control characters.

## Render lineages

Renders the complete lineage of a term

    bio taxon 564286 --lineage

prints:

    no rank, root, 1
       no rank, cellular organisms, 131567
          superkingdom, Bacteria (eubacteria), 2
             clade, Terrabacteria group, 1783272
                phylum, Firmicutes, 1239
                   class, Bacilli, 91061
                      order, Bacillales, 1385
                         family, Bacillaceae, 186817
                            genus, Bacillus, 1386
                               species group, Bacillus subtilis group, 653685
                                  species, Bacillus subtilis, 1423
                                     strain, Bacillus subtilis str. 10, 564286

## Filter blast results

(TODO) - filters BLAST alignment for species that fall within a taxonomical clade

## Preloading data

In many usecases,  the default behavior is extremely fast and can produce family, genus and species level information in fractions of a second.

But when the number of nodes is very large, while the output starts streaming immediately the software will query the database for each child node separately, hence making the process slower:

For example the command below will render the complete NCBI tree with over 2 million entries and might take a very long time to produce (15 minutes or so)

    bio taxon 1

The system offers a speedup process of preloading all the data into memory to speed up the process at the cost of imposing a 6 second loading penalty.

    bio taxon 1 --preload
    
The command above takes a total of just 10 seconds to generate a large tree with the entire NCBI taxonomical tree.    

For queries that take more than 10 seconds to complete we recommend applying the `--preload` flag.

 