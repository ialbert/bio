# Fetch sequence from NCBI
bio fetch NC_045512 MN996532 > genomes.gb

# Convert the sequences of interest
bio convert genomes.gb --seqid Wuhan:S --fasta  > s1.fa

