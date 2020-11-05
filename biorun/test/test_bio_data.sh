#
# Build the new test sets when the output formats change.
#

# Stop on errors.
set -uex

# Delete the SARS2 data if exists
bio SARS2 --delete

# Rebuild the JSON
bio NC_045512 --fetch --rename SARS2 --seqid SARS2

# JSON output.
bio SARS2 > SARS2.json

# FASTA.
bio SARS2 --fasta > SARS2.fa

# GFF formatting.
bio SARS2 --gff > SARS2.gff

# JSON by type and match
bio SARS2 --match ORF1ab --type gene > match.json

# GFF by gene name.
bio SARS2 --gff --gene S > gene.gff

# GFF by start and end.
bio SARS2 --gff  --start 10000 --end 20000 > overlap.gff

# GFF by type.
bio SARS2 --gff  --type CDS > type.gff

# Sliced FASTA with different id .
bio SARS2 --fasta --seqid foo --start 10 --end 20 > fasta-start.fa

# FASTA features
bio SARS2 --fasta --type CDS > CDS.fa

# Sliced FASTA features by type.
bio SARS2 --fasta --type gene --end 10 > gene-start.fa

# Protein FASTA.
bio SARS2 --protein --start -10 > protein-end.fa

# Renamed protein.
bio SARS2:S --fasta --protein --seqid foo > s_prot_foo.fa
