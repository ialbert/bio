set -uex

ACC=NC_045512

GBK=$ACC.gb

XML=$ACC.xml

# Store the output here.
mkdir -p out

# Get the input file as GenBank
#(cd out && efetch -db nuccore -id $ACC -format gbwithparts > $ACC.gb)

# Get the input file as FASTA
#(cd out && efetch -db nuccore -id $ACC -format fasta > $ACC.fa)

# Get the input file as XML
#(cd out && efetch -db nuccore -id $ACC -format gbwithparts -mode xml > $ACC.xml)

# Seqret conversion
#(cd out && cat $GBK | seqret -filter -feature -osformat gff3 > $ACC.seqret.gff)

# Readseq conversion
#(cd out && cat $GBK | readseq -p -format=GFF  > $ACC.readseq.gff)

# Bioperl conversion
#(cd out && cat $GBK | bp_genbank2gff3.pl -in stdin -out stdout > $ACC.bioperl.gff)

# Biocode bassed conversion
#(cd out && convert_genbank_to_gff3.py -i $GBK -o $ACC.biocodes.gff)

# BCCB Biopython conversion
#(cd out && python ../genbank_to_gff.py $GBK > $ACC.bcbb.gff)

# Lindenbaum's conversion (NULL Pointer Exception!)
#(cd out && cat $XML | java -jar ../gb2gff.jar > $ACC.varkit.gff)

# Use the bio conversion
bio $ACC --fetch --gff

