# bio format: summarize alignments {#bio-format}

Note: the inputs for this tool must be so in the so called "aligned" FASTA. The sequences need to arranged to contain gaps in the places of insertions and deletions. For example:


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

    mafft --auto --preservecase genomes.fa > aligned.fa

Generate the variant file:

    cat aligned.fa | bio format > variants.vcf

investigate the variants:

    cat variants.vcf | head

prints:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">
    ##contig=<ID=MN996532.2,length=29875,assembly=MN996532.2>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NC_045512.2
    MN996532.2	19	19_T_C	T	C	.	PASS	TYPE=SNP	GT	1
    MN996532.2	37	37_G_C	G	C	.	PASS	TYPE=SNP	GT	1
    MN996532.2	42	42_C_T	C	T	.	PASS	TYPE=SNP	GT	1
    MN996532.2	91	91_A_G	A	G	.	PASS	TYPE=SNP	GT	1


## One shot wonder

You can do the above in a single command:

    bio fetch NC_045512 MN996532 | bio fasta | mafft --preservecase - | bio format > variants.vcf



