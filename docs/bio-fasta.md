# GenBank to FASTA {#bio-fasta}

A GenBank file represents sequence information in multiple sections:

1. Genomic sequences (the entire genomic sequence)
1. Feature annotation (intervals relative to the genome)

    
### Get a GenBank file

    bio fetch NC_045512,MN996532 > genomes.gb

### Convert to FASTA

The default behavior is to convert the genome (source) of the GenBank file to FASTA. The following commands print the genome sequence:

    bio fasta genomes.gb

or it also works as a stream
    
    cat genomes.gb | bio fasta

### Selecting features

If any feature selector is passed the FASTA conversion operates on the features in the GenBank:

    bio fasta --type CDS genomes.gb

will convert to fasta the coding sequences alone.


### Get the sequence for the genome

    bio fasta genomes.gb | head -3

### Manipulate a genomic subsequence

    bio fasta genomes.gb --start 100 --end 10kb 

### Extract the sequences for annotations of a certain type

    bio fasta genomes.gb --type CDS | head -3

### Extract CDS sequences by gene name

    bio fasta --gene S--end 60 genomes.gb 

### Extract sequence by feature accession number

    cat genomes.gb | bio fasta --end 10 --name QHR63308.1

### Translate the sequence

This command translates the DNA sequence to peptides:

    bio fasta --end 30 --translate --gene S genomes.gb

The slice to 30 is applied on the DNA sequence before the translation.

### Extract the protein sequence

This flag extracts the protein sequence embedded in the original GenBank file:

    bio fasta --end 10 --protein --gene S genomes.gb

Note how in this case the slice to 10 is applied on the protein sequence.
