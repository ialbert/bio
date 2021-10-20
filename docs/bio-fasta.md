# bio fasta: convert to FASTA {#bio-fasta}

Install `bio` with:

    pip install bio --upgrade

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.

## Rationale

GenBank/EMBL files represents sequence information in multiple sections:

1. Genomic sequences (the entire genomic sequence)
1. Feature annotation (intervals relative to the genome)

`bio fasta` can operate on GenBank/EMBL files, filter and extract various subsets of the data.

## Get a GenBank file

    bio fetch NC_045512 MN996532 > genomes.gb

## Convert to FASTA

The default behavior is to convert the genome the GenBank file to FASTA:

    bio fasta genomes.gb > genomes.fa

to convert the features  component pass the `--features` flag or use any of `--type`, `--gene` or other feature specific selectors.

    bio fasta genomes.gb --features > features.fa

## Inputs

The input may be GENBANK, FASTA, EMBL or FASTQ.

## What gets converted?

GenBank and EMBL files contain both genomes and features all features are extracted.

    cat genomes.gb | bio fasta > genomes.fa

pass any feature matcher to limit to certain types:

    bio fasta --type CDS genomes.gb -e 10 | head

prints:

    >YP_009724389.1 {"type": "CDS", "gene": "ORF1ab", "product": "ORF1ab polyprotein", "locus": "GU280_gp01"}
    ATGGAGAGCC
    >YP_009725295.1 {"type": "CDS", "gene": "ORF1ab", "product": "ORF1a polyprotein", "locus": "GU280_gp01"}
    ATGGAGAGCC

## Select by name

`-m` or `--match` performs a regular expression match on sequence ids:

    cat genomes.gb | bio fasta -m glyco -end 10

prints:

    >YP_009724390.1 {"type": "CDS", "gene": "S", "product": "surface glycoprotein", "locus": "GU280_gp02"}
    ATGTTTGTTT
    >YP_009724393.1 {"type": "CDS", "gene": "M", "product": "membrane glycoprotein", "locus": "GU280_gp05"}T

`-i` or `--id` performs an exact match on sequence ids:

    cat genomes.gb | bio fasta -i YP_009724390.1  -end 10

prints:

    >YP_009724390.1 {"type": "CDS", "gene": "S", "product": "surface glycoprotein", "locus": "GU280_gp02"}
    ATGTTTGTTT

pass multiple ids to match multiple sequences:

    cat genomes.gb | bio fasta -i YP_009724390.1,QHR63300.2  -end 10

prints:

    >YP_009724390.1 {"type": "CDS", "gene": "S", "product": "surface glycoprotein", "locus": "GU280_gp02"}
    ATGTTTGTTT
    >QHR63300.2 {"type": "CDS", "gene": "S", "product": "spike glycoprotein", "locus": ""}
    ATGTTTGTTT

## Selecting features

If any feature selector is passed the FASTA conversion operates on the features in the GenBank:

    bio fasta --type CDS genomes.gb

will convert to fasta the coding sequences alone.


## Manipulate a genomic subsequence

    bio fasta genomes.gb --start 100 --end 10kb

## Extract the sequences for annotations of a certain type

    bio fasta genomes.gb --type CDS | head -3

### Extract CDS sequences by gene name

    bio fasta --gene S--end 60 genomes.gb 

### Extract sequence by feature accession number

    cat genomes.gb | bio fasta --end 10 --id QHR63308.1

### Translate the sequence

This command translates the DNA sequence to peptides:

    bio fasta --end 30 --translate --gene S genomes.gb

The slice to 30 is applied on the DNA sequence before the translation.

### Extract the protein sequence

This flag extracts the protein sequence embedded in the original GenBank file:

    bio fasta --end 10 --protein --gene S genomes.gb

Note how in this case the slice to 10 is applied on the protein sequence.

## Other tools 

My first choice when needing functionality not present in `bio fasta` would be to look at:

* `seqkit` a cross-platform and ultrafast toolkit for FASTA/Q file manipulation <https://bioinf.shenwei.me/seqkit/>

and

* `seqtk` a fast and lightweight tool for processing sequences in the FASTA or FASTQ format. It seamlessly parses both FASTA and FASTQ files which can also be optionally compressed by gzip.  <https://github.com/lh3/seqtk>


### Other potentially useful software

The following software may be installed with `conda/mamba`:

* `any2fasta` converts various outputs to FASTA see <https://github.com/tseemann/any2fasta>
* `pyfasta` pythonic access to fasta sequence files <https://github.com/brentp/pyfasta>
* `pyfastq` Python3 script to manipulate FASTA and FASTQ (and other format) files, plus API for developers https://github.com/sanger-pathogens/Fastaq   
* `seqret` part the EMBOSS suite <https://www.bioinformatics.nl/cgi-bin/emboss/help/seqret>
* `bam2fasta` Convert 10x bam file or fastq.gz files to individual FASTA files per cell barcode convert large bam files to fastq.gz format before the individual fasta files per cell barcode conversion <https://github.com/czbiohub/bam2fasta>
