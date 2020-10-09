#
# Build the new test sets when the output formats change.
#

# Stop on errors.
set -uex

# Complete GFF formatting.
bio view NC_045512 --fasta > NC_045512.fa

# Complete GFF formatting.
bio view NC_045512 --gff > NC_045512.gff

# Selection by gene name.
bio view NC_045512 --gff --gene S > parts/gene.gff

# Selection by start and end.
bio view NC_045512 --gff  --start 21563 --end 21565 > parts/start.gff

# Selection by type.
bio view NC_045512 --gff  --type CDS > parts/type.gff