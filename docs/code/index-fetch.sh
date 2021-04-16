# Fetch sequence from NCBI
bio fetch NC_045512 MN996532 > genomes.gb

# Covert the sequences into a JSON format
cat genomes.gb | bio json > genomes.json

# You can also do it in one step
bio fetch NC_045512 MN996532 | bio json > genomes.gb
