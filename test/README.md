Demo accession numbers:

Viral, SARS-COV-2: NC_045512 (1 second)

    time bio NC_045512 --fetch

Viral, Ebola: AF086833 (1 second)

    time bio AF086833  --fetch

Another Ebola (2014):

    time bio KJ660346 --fetch --rename ebola14

# Transcript 
/transcript_id="NM_164349.3"
/db_xref="FLYBASE:FBtr0078171"
/db_xref="GeneID:33156"
/db_xref="FLYBASE:FBgn0002121"

# Protein
/protein_id="NP_722583.1"
/db_xref="FLYBASE:FBpp0077824"
/db_xref="GeneID:33156"
/db_xref="FLYBASE:FBgn0002121"

Bacteria, Escherichia coli: NC_002695 (2 minutes, 10MB)

    time bio NC_002695 --fetch -v 

Insect, fruit fly (chromosome 2L): NT_033779 (2 minutes, 44MB)

    time bio NT_033779 --fetch 

Mammal, human genome (chromosome 1): NC_000001 (2 minutes, 300MB)

    time bio NC_000001  --fetch
    

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
    

   