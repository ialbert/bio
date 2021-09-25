set -uex

ACC=NC_045512
VER=2

GBK=$ACC.gb

XML=$ACC.xml

# The official SARS-COV-2 annotation
#URL='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz'


# Store the output here.
mkdir -p out

# Get the input file as GenBank
#(cd out && efetch -db nuccore -id $ACC -format gbwithparts > $ACC.gb)

# Get the input file as FASTA
#(cd out && efetch -db nuccore -id $ACC -format fasta > $ACC.fa)

# Get the input file as XML
#(cd out && efetch -db nuccore -id $ACC -format gbwithparts -mode xml > $ACC.xml)

#(cd out && curl $URL | gunzip -c > $ACC.ncbi.gff)

# Seqret conversion
#(cd out && cat $GBK | seqret -filter -feature -osformat gff3 > $ACC.seqret.gff)

# Readseq conversion
#(cd out && cat $GBK | readseq -p -format=GFF  > $ACC.readseq.gff)

# Bioperl conversion. Adds spurious first line!
(cd out && bp_genbank2gff3.pl $GBK -out stdout | python ../fixit.py $ACC $VER > $ACC.bioperl.gff)

# Biocode bassed conversion
#(cd out && convert_genbank_to_gff3.py -i $GBK -o $ACC.biocodes.gff)

# BCCB Biopython conversion, we need to rename the sequence to the locus.
#(cd out && python ../genbank_to_gff.py $GBK | awk -v ACC=$ACC ' !/#/ { $1=ACC; print $0 }' > $ACC.bcbb.gff)

# Lindenbaum's conversion (NULL Pointer Exception!)
#(cd out && cat $XML | java -jar ../gb2gff.jar > $ACC.varkit.gff)

# Use the bio conversion
(cd out && bio $ACC --fetch --gff > $ACC.bio.gff)

