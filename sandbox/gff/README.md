A collection of Genbank to GFF conversion scripts

Naively I was hoping that the problem of converting GenBank to GFF has been solved ... I mean after all these years ... turns out it is not at all.

How many ways are there to  to convert GenBank to GFF? Too many. Below is what I found in one afternoon.

After trying out and benchmarking each turns out each has serious problems that make them unusable.

See the `runall.sh` script

# Bioperl
https://github.com/bioperl/bioperl-live/blob/master/bin/bp_genbank2gff3

# varkit : Java utilities for Bioinformatics
http://lindenb.github.io/jvarkit/GenbankToGff3.html

# BCBB, Biopython 
https://github.com/chapmanb/bcbb/tree/master/gff

# Biocode based
https://github.com/jorvis/biocode
https://github.com/jorvis/biocode/blob/master/gff/convert_genbank_to_gff3.py

pip install biocode

# EMBOSS
http://emboss.sourceforge.net/apps/cvs/emboss/apps/seqret.html

# ReadSeq
cat NC_001501.gb | readseq -p -format=FASTA > NC_001501-version3.fa
