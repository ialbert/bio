# Convert to FASTA {#bio-fasta}

A GenBank file represents sequence inforamtion in multiple ways:

1. Genomic sequences (the entire genomic sequence)
1. Feature annotation (intervals relative to the genome)

In `bio` we operate on:

* `--genome` to access the genomic sequence
* `--fasta` to access the feature sequences

both output are in Fasta format.

### Get a dataset

Get SARS-COV-2 data and rename it to `ncov`:

```{bash, comment=NA}
bio NC_045512 --fetch --rename ncov
```

### Get the sequence for the genome

```{bash, comment=NA}
bio ncov --genome | head -3
```

### Manipulate a genomic subsequence

```{bash, comment=NA}
bio ncov --genome --start 100 --end 130 --seqid foo 
```


### Extract the sequences for annotations of a certain type

```{bash, comment=NA}
bio ncov --fasta --type CDS | head -3
```

### Extract CDS sequences by gene name

```{bash, comment=NA}
bio ncov --gene S --fasta --end 60 
```

a shortcut notation of the above:

```{bash, comment=NA}
bio ncov:S --fasta --start 100 --end 150 
```

### Extract sequence by feature accession number

```{bash, comment=NA}
bio ncov -id YP_009724390.1 --fasta --start 100 --end 150 
```

### Translate the sequence

This command translates the DNA sequence to peptides:

```{bash, comment=NA}
bio ncov:S --fasta --end 180 --translate
```

The slice to 180 is applied on the DNA sequence before the translation.

### Extract the protein sequence

This flag extracts the protein sequence embedded in the original GenBank file:

```{bash, comment=NA}
bio ncov:S --fasta --end 60 --protein
```

Note how in this case the slice to 60 is applied on the protein sequence.