#
# Build the new test sets when the output formats change.
#

# Stop on errors.
set -uex

# Delete the ncov data if exists
bio ncov --delete

# Rebuild the JSON
bio NC_045512 --fetch --rename ncov --seqid ncov

# JSON output.
bio ncov > ncov.json

# FASTA.
bio ncov --fasta > ncov.fa

# GFF formatting.
bio ncov --gff > ncov.gff

# JSON by type and match
bio ncov --match ORF1ab --type gene > match.json

# GFF by gene name.
bio ncov --gff --gene S > gene1.gff

# GFF by start and end.
bio ncov --gff  --start 10000 --end 20000 > overlap.gff

# GFF by type.
bio ncov --gff  --type CDS > type.gff

# Sliced FASTA with different id .
bio ncov --fasta --seqid foo --start 10 --end 20 > fasta-start.fa

# FASTA features
bio ncov --fasta --type CDS > CDS.fa

# Sliced FASTA features by type.
bio ncov --fasta --type gene --end 10 > gene-start.fa

# Protein FASTA.
bio ncov --protein --start -10 > protein-end.fa

# Translated CDS.
bio ncov --translate --type CDS > translate.fa

# Renamed protein.
bio ncov:S --fasta --protein --seqid foo > s_prot_foo.fa

# Get the RatG13 data.
bio MN996532 --fetch --rename ratg13 --seqid ratg13

# Align DNA
bio align ncov ratg13 --end 200 > align-dna.txt

# Align the extracted protein.
bio align ncov:S ratg13:S --end 80 > align-dna-s.txt

# Align the extracted protein.
bio align ncov:S ratg13:S --protein > align-protein-s.txt

# Align the translated regions.
bio align ncov:S ratg13:S --end 80 --translate > align-translated-s.txt
