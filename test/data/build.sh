#
# Build the new test sets when the output formats change.
#

# Stop on errors.
set -uex

# Rebuild the JSON
bio fetch NC_045512 --build  --name SARS2

# JSON output.
bio view SARS2 > SARS2.json

# JSON by type and match
bio view SARS2 --match ORF1ab --type gene > parts/match.json

# GFF formatting.
bio view SARS2 --gff > SARS2.gff

# GFF by gene name.
bio view SARS2 --gff --gene S > parts/gene.gff

# GFF by start and end.
bio view SARS2 --gff  --start 10000 --end 20000 > parts/overlap.gff

# GFF by type.
bio view SARS2 --gff  --type CDS > parts/type.gff

# FASTA.
bio view SARS2 --fasta > SARS2.fa

# Sliced FASTA with different id .
bio view SARS2 --fasta --id foo --start 10 --end 20 > parts/fasta-start.fa

# FASTA features
bio view SARS2 --fasta --type CDS > parts/CDS.fa

# Sliced FASTA features by type.
bio view SARS2 --fasta --type gene --end 10 > parts/gene-start.fa

# Protein FASTA
bio view SARS2 --protein --start -10 > parts/protein-end.fa
