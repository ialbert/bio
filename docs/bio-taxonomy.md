# bio taxon: show taxonomies {#bio-taxonomy}

The `bio` package provides utilities to visualize NCBI taxonomies. See the help with:

    bio taxon -h

## Download or build the taxonomy

Before using the taxonomy related database needs to downloaded with:

    bio --download 


## Check database

    bio taxon
    
prints:

    Content: 2,359,690 taxon names; 209,096 parent ranks

There are a total of `2,353,384 ` nodes (taxonomical entries) out of which `207,447` are nodes that are non-terminal (non-leaf) nodes. For these numbers we see that the vast majority of the taxonomy annotations are for terminal, leaf nodes.

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

The search words may use regular expression control characters, selecting species names starting with `jawed`:

    bio taxon '^jawed'

produces:

    # searching taxonomy for: ^jawed
    clade, Gnathostomata (jawed vertebrates), 7776

## View taxonomy for data 

Once you fetch the data
    
    bio fetch NC_045512 > genome.gb
        
you can view the lineage of the taxonomical id in the GenBank

    bio taxon genome.gb --lineage

that prints:

    superkingdom, 10239, Viruses
      clade, 2559587, Riboviria
        kingdom, 2732396, Orthornavirae
          phylum, 2732408, Pisuviricota
            class, 2732506, Pisoniviricetes
              order, 76804, Nidovirales
                suborder, 2499399, Cornidovirineae
                  family, 11118, Coronaviridae
                    subfamily, 2501931, Orthocoronavirinae
                      genus, 694002, Betacoronavirus
                        subgenus, 2509511, Sarbecovirus
                          species, 694009, Severe acute respiratory syndrome-related coronavirus
                            no rank, 2697049, Severe acute respiratory syndrome coronavirus 2

## View taxonomy by tax id
    
Pass a NCBI taxonomical id to see all the descendants of it:

    bio taxon 117565 | head

prints:

    class, 117565, Myxini
      order, 7761, Myxiniformes
        family, 7762, Myxinidae
          subfamily, 30309, Eptatretinae
            genus, 7763, Eptatretus
              species, 7764, Eptatretus burgeri
              species, 7765, Eptatretus stoutii
              species, 7767, Eptatretus okinoseanus
              species, 50612, Eptatretus atami
              species, 78394, Eptatretus cirrhatus

we can limit the depth (how many levels to show):

    bio taxon 117565 -d 3

prints:

    class, 117565, Myxini
      order, 7761, Myxiniformes
        family, 7762, Myxinidae
        no rank, 727487, unclassified Myxiniformes

## Print lineage

    bio taxon 564286 --lineage

prints:

    no rank, 131567, cellular organisms
      superkingdom, 2, Bacteria
        clade, 1783272, Terrabacteria group
          phylum, 1239, Firmicutes
            class, 91061, Bacilli
              order, 1385, Bacillales
                family, 186817, Bacillaceae
                  genus, 1386, Bacillus
                    species group, 653685, Bacillus subtilis group
                      species, 1423, Bacillus subtilis
                        strain, 564286, Bacillus subtilis str. 10

## Filter results by taxonomic clades

In bioinformatics we frequently have files where a column corresponds to taxonomic ids and we wish to filter rows
by clades. For example keep all rows where the taxid is belongs to a subtree in the taxonomy. Here is a file

    11138,AB008939.1
    11138,AB008940.1
    215681,AB242262.1
    215681,AB242263.1
    215681,AB242264.1
    229992,AB257344.1
    391355,AB263618.1
    11128,AB277098.1
    422141,AB277099.1
    11133,AB277100.1

We can summarize the taxonomies in the file with:

    cat data.txt | bio taxon -d 1

to print:

    no rank,11138,Murine hepatitis virus
    no rank,215681,Canine respiratory coronavirus
    no rank,229992,SARS coronavirus Frankfurt 1
    no rank,391355,SARS coronavirus Frankfurt1-v01
    no rank,11128,Bovine coronavirus
    no rank,422141,Giraffe coronavirus US/OH3-TC/2006
    no rank,11133,Bovine coronavirus strain Quebec

Which data are a kind of "Bovine Coronavirus"?

    cat data.txt | bio taxon -d 1 -keep 11128

prints:

    11128,AB277098.1
    422141,AB277099.1
    11133,AB277100.1

Which data are not a kind of "Bovine Coronavirus"?

    cat data.txt | bio taxon -d 1 -remove 11128

prints:

    11138,AB008939.1
    11138,AB008940.1
    215681,AB242262.1
    215681,AB242263.1
    215681,AB242264.1
    229992,AB257344.1
    391355,AB263618.1

## Chain the actions

    cat data.txt | bio taxon -d 1 -remove 11128 | bio taxon -d 1

prints:

    no rank,11138,Murine hepatitis virus
    no rank,215681,Canine respiratory coronavirus
    no rank,229992,SARS coronavirus Frankfurt 1
    no rank,391355,SARS coronavirus Frankfurt1-v01

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

Note: this command benefits greatly from using `--preload`. This command will only finish in reasonable time (in 10 seconds) with preload. Since the entire database is visited, it has to be loaded into memory for it to be fast.

    bio taxon --list --preload | wc -l

## Update the taxonomy

You may build the newest version locally:

    bio taxon --build
    
The command will download and build a new taxonomy using the latest NCBI taxonomy data. The efficiency of the process depends on the speed of the hard drive and takes around 20 minutes.

## Preloading data

For many use cases,  the default behavior is plenty fast and can produce family, genus and species level information in a fraction of a second.

Internally, during operation, the software will query the database for each child node. When selecting a rank where the number of descendant nodes is large (over 10,000 nodes) the run time of the independent queries adds up to a substantial overhead.

For example the command below attempts to render the complete NCBI taxonomic tree with over 2.2 million (!) descendant nodes. When run like so it will take a very long time to produce the output (more than two hours):

    bio taxon 1  | head

The data starts scrolling immediately, but the individual queries add up to potentially hours of runtime.

The software can operate in a different mode to speed up the process by preloading all the data into memory at the cost of imposing a 8 second pre-loading penalty.

    bio taxon 1 --preload | wc -l
    
When run with the `--preload` flag the command takes a total of just 11 seconds to generate the same large tree of the entire NCBI taxonomical tree. We don't apply this mode by default because all queries would then take at least 6 seconds, even those that currently finish very quickly.

For queries that take more than 10 seconds to complete (have more than 10,000 descendant nodes) we recommend applying the `--preload` flag.

## Build the newest version locally

You may build the newest version locally, will take about an hour:

    bio taxon --build
