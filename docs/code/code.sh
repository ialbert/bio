set -uex

bio align GATTACA GATCA  > align1.txt

bio align GATTACA GATCA --table | column -t > align2.txt

bio align GATTACA GATCA --vcf > align3.txt

bio align GATTACA GATCA --diff | column -t > align4.txt

bio align GATTACA GATCA --global > align5.txt

#bio fetch MN996532 NC_045512 > genomes.gb

cat genomes.gb | bio fasta --gene S | bio align | head -8 > align6.txt

cat genomes.gb | bio fasta --gene S --protein | bio align --diff | head -8 > align7.txt

cat genomes.gb | bio fasta --gene S  --protein | bio align --table | column -t > align8.txt

cat genomes.gb | bio fasta --gene S --protein | bio align --diff | tail -5 | column -t > align9.txt

cat genomes.gb | bio fasta --gene S --translate | bio align --matrix PAM30 | head -8 > align10.txt

cat genomes.gb | bio fasta > genomes.fa

# Alignments

mafft --auto --quiet --preservecase genomes.fa  > aligned.fa

cat aligned.fa | bio format | head -12  > format1.txt

cat aligned.fa | bio format --diff | head -5 > format2.txt

cat aligned.fa | bio format --vcf | head > format3.txt

#
# Table
#
cat genomes.gb | bio table | head -5 | column -t > table1.txt

cat genomes.gb | bio table --fields id,type,size,gene | head -5 | column -t > table2.txt

cat genomes.gb | bio table --type CDS --fields id,gene,isolate,country,date | head -5 | column -t > table3.txt

#
# Search examples
#
bio search AF086833 > search1.txt

