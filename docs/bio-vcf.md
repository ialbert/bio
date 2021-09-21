# bio vcf: VCF from aligned FASTA {#bio-vcf}

Note: the inputs for this tool must be so called "aligned" FASTA. The sequences need to arranged to contain gaps in the place of insertions and deletions. For example:


    >seq1
    THIS-LI-NE-
    >seq2
    --ISALIGNED

Many tools can generate such alignments.

    mamba install mafft


## Usage example

Get data for coronavirus and a related bat coronavirus

    bio fetch NC_045512 MN996532 > genomes.gb

Turn the genomes into FASTA

    cat genomes.gb | bio fasta > genomes.fa


Align the fasta files with `mafft`:

    mafft --auto genomes.fa > aligned.fa


Generate the variant file:

    cat aligned.fa | bio fasta2vcf > variants.vcf

investigate the variants:

    cat variants.vcf | head

prints:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">
    ##contig=<ID=NC_045512.2,length=29903,assembly=NC_045512.2>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MN996532.2
    NC_045512.2	19	C19T	C	T	.	PASS	TYPE=SNP	GT	1
    NC_045512.2	37	C37G	C	G	.	PASS	TYPE=SNP	GT	1
    NC_045512.2	42	T42C	T	C	.	PASS	TYPE=SNP	GT	1
    NC_045512.2	91	G91A	G	A	.	PASS	TYPE=SNP	GT	1


## One shot wonder

You can do the above in a single command:

    bio fetch NC_045512 MN996532 | bio fasta | mafft - | bio fasta2vcf > variants.vcf



