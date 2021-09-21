Demo accession numbers:

Viral, SARS-COV-2: NC_045512 (1 second)

    time bio fetch NC_045512

Viral, Ebola: AF086833 (1972), KJ660346 (2014)

    time bio fetch AF086833  KJ660346  > genome.gb

# Transcript 


Bacteria, Escherichia coli: NC_002695 (2 minutes, 10MB)

    time bio fetch  NC_002695

Insect, fruit fly (chromosome 2L): NT_033779 (2 minutes, 44MB)

    time bio fetch NT_033779

Mammal, human genome (chromosome 1): NC_000001 (2 minutes, 300MB)

    time bio fetch NC_000001

A weird join operator:

    NC_003977


Multiple equally scoring alignments:

    bio align AGGATAGGAG AGGATTAG -i -v

Translating global alignemnt:

    bio align AATATTAGA AGATGAGAG -i -pep1
 
# Taxonomy test

    # Flowering plant
    bio taxon 70730 --lineage
     
    # Bacillus subtilis 
    bio taxon 564286 --lineage
    
    # Yucatan goby
    bio taxon 70730 --lineage
         
Complicated GO term:

    bio --define intergenic mrna trans splicing
    

    # Multiple equally scoring alignments
    bio fasta genomes.gb --gene S -s 1432 --end 1500 | bio align


