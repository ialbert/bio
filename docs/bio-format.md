# bio format: summarize alignments {#bio-format}

Note: the inputs for this tool must be in the so called "aligned" FASTA. The sequences need to arranged to contain gaps in the places of insertions and deletions. For example:


    >seq1
    THIS-LI-NE-
    >seq2
    --ISALIGNED

Many tools can generate such alignments.

    mamba install mafft

## Usage example

Get FASTA sequences for the coronavirus and the closest related bat coronavirus

    bio fetch NC_045512 MN996532 | bio fasta --genome > genomes.fa

Align the fasta files with `mafft`:

    mafft --auto --quiet --preservecase genomes.fa  > aligned.fa

## Pairwise alignment

The default output generates the diffs:

    cat aligned.fa | bio format | head -12

shows the pairwise alignment:

    # NC_045512.2 (29903) vs MN996532.2 (29855)
    # pident=96.0% len=29903 ident=28720 mis=1135 del=0 ins=48

    ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAA
    ||||||||||||||||||.|||||||||||||||||.||||.|||||||||||||||||||||||||||||||||||||||
    ATTAAAGGTTTATACCTTTCCAGGTAACAAACCAACGAACTCTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAA

    ATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGAC
    |||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ATCTGTGTGACTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGAC


## Show differences

    cat aligned.fa | bio format --diff | head -5

This output generates an easy-to-read output of the alignments:

    T19C    SNP     19      T       C
    G37C    SNP     37      G       C
    C42T    SNP     42      C       T
    A91G    SNP     91      A       G
    A174G   SNP     174     A       G

Columns are: position, type, query, change, target name.

## Generate VCF output

    cat aligned.fa | bio format --vcf | head

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

    bio fetch NC_045512 MN996532 | bio fasta | mafft --preservecase - | bio format --vcf > variants.vcf



