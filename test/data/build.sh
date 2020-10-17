#
# Build the new test sets when the output formats change.
#

# Stop on errors.
set -uex

# Rebuild the JSON
bio fetch NC_045512 --build

# JSON output.
bio view NC_045512 > NC_045512.json

# JSON by type and match
bio view NC_045512 --match ORF1ab --type gene > parts/match.json

# GFF formatting.
bio view NC_045512 --gff > NC_045512.gff

# GFF by gene name.
bio view NC_045512 --gff --gene S > parts/gene.gff

# GFF by start and end.
bio view NC_045512 --gff  --start 10000 --end 20000 > parts/overlap.gff

# GFF by type.
bio view NC_045512 --gff  --type CDS > parts/type.gff

# FASTA.
bio view NC_045512 --fasta > NC_045512.fa

# Sliced FASTA with different id .
bio view NC_045512 --fasta --id foo --start 10 --end 20 > parts/fasta-start.fa

# FASTA features
bio view NC_045512 --fasta --type CDS > parts/CDS.fa

# Sliced FASTA features by type.
bio view NC_045512 --fasta --type gene --end 10 > parts/gene-start.fa

# Protein FASTA
bio view NC_045512 --protein --start -10 > parts/protein-end.fa
