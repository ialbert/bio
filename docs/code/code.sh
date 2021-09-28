
bio align GATTACA GATCA  > align1.txt

bio align GATTACA GATCA --table | column -t > align2.txt


bio align GATTACA GATCA --vcf > align3.txt

bio align GATTACA GATCA --diff | column -t > align4.txt

bio align GATTACA GATCA --global > align5.txt

bio fetch NC_045512 MN996532 > genomes.gb

cat genomes.gb | bio fasta --gene S | bio align | head -8> align6.txt

cat genomes.gb | bio fasta --gene S --translate | bio align | head -8 > align7.txt

cat genomes.gb | bio fasta --gene S --translate | bio align --table | column -t > align8.txt

cat genomes.gb | bio fasta --gene S --translate | bio align --diff | column -t | tail -5 > align9.txt

cat genomes.gb | bio fasta --gene S --translate | bio align --matrix PAM30 | head -8 > align10.txt
