set -uex

bio align GATTACA GATCA  > align1.txt

bio align GATTACA GATCA --table | column -t > align2.txt

bio align GATTACA GATCA --vcf > align3.txt

bio align GATTACA GATCA --mut | column -t > align4.txt

bio align GATTACA GATCA --global > align5.txt

#bio fetch MN996532 NC_045512 > genomes.gb

cat genomes.gb | bio fasta --gene S | bio align | head -8> align6.txt

cat genomes.gb | bio fasta --gene S --protein | bio align --mut | head -8 > align7.txt

cat genomes.gb | bio fasta --gene S --translate | bio align --table | column -t > align8.txt

cat genomes.gb | bio fasta --gene S --translate | bio align --matrix PAM30 | head -8 > align10.txt

cat genomes.gb | bio fasta --genome > genomes.fa

mafft --auto --quiet --preservecase genomes.fa  > aligned.fa

cat aligned.fa | bio format | head -12  > format1.txt

cat aligned.fa | bio format --mut | head -5 > format2.txt

cat aligned.fa | bio format --vcf | head > format3.txt
