# Convert to FASTA {#bio-fasta}

A GenBank file represents sequence inforamtion in multiple ways:

1. Genomic sequences (the entire genomic sequence)
1. Feature annotation (intervals relative to the genome)

In `bio` we operate on:

* `--fasta` to access the genome
* `--fasta --features` to access the features annotated on the genome

The `--features` flag often  not necessary as   `bio` will set it automatically if it is obvious that the command targets features. For example `--type CDS` will turn on feature rendering mode.

## Shortcuts

A a colon delimited term:

    bio convert foo:bar --fasta

is equivalent to: 

    bio convert foo --type CDS --gene bar --fasta
    
### Get a dataset

Get SARS-COV-2 data and rename it to `ncov`:

```{bash, comment=NA}
bio fetch NC_045512 --rename ncov
```

### Get the sequence for the genome

```{bash, comment=NA}
bio convert ncov --fasta | head -3
```

### Manipulate a genomic subsequence

```{bash, comment=NA}
bio convert ncov --fasta --start 100 --end 130 --seqid foo 
```


### Extract the sequences for annotations of a certain type

```{bash, comment=NA}
bio convert ncov --fasta --type CDS | head -3
```

### Extract CDS sequences by gene name

```{bash, comment=NA}
bio convert ncov --gene S --fasta --end 60 
```

a shortcut notation of the above:

```{bash, comment=NA}
bio convert ncov:S --fasta --start 100 --end 150 
```

### Extract sequence by feature accession number

```{bash, comment=NA}
bio convert ncov -id YP_009724390.1 --fasta --start 100 --end 150 
```

### Translate the sequence

This command translates the DNA sequence to peptides:

```{bash, comment=NA}
bio convert ncov:S --fasta --end 180 --translate
```

The slice to 180 is applied on the DNA sequence before the translation.

### Extract the protein sequence

This flag extracts the protein sequence embedded in the original GenBank file:

```{bash, comment=NA}
bio convert ncov:S --fasta --end 60 --protein
```

Note how in this case the slice to 60 is applied on the protein sequence.