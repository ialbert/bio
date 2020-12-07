Demo accession numbers:

Viral, SARS-COV-2: NC_045512 (1 second)

    time bio NC_045512 --fetch

Viral, Ebola: AF086833 (1 second)

    time bio AF086833  --fetch

Another Ebola (2014):

    time bio KJ660346 --fetch --rename ebola14

Bacteria, Escherichia coli: NC_002695 (2 minutes, 10MB)

    time bio NC_002695 --fetch -v 

Insect, fruit fly (chromosome 2L): NT_033779 (2 minutes, 44MB)

    time bio NT_033779 --fetch 

Mammal, human genome (chromosome 1): NC_000001 (2 minutes, 300MB)

    time bio NC_000001  --fetch
    

A weird join operator:

    NC_003977


Multiple equally scoring alignments:

    bio AGGATAGGAG AGGATTAG -i --align -v

# Taxonomy test

    # Flowering plant
    bio taxon 70730 --lineage
     
    # Bacillus subtilis 
    bio taxon 564286 --lineage
    
    # Yucatan goby
    bio taxon 70730 --lineage
         