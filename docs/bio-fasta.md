# Convert to FASTA {#bio-fasta}

### Get a dataset

Get SARS-COV-2 data and rename it to `ncov`:

```{bash, comment=NA}
    bio NC_045512 --fetch --rename ncov
```

## Convert the data to FASTA (the origin of the GenBank)

```{bash, comment=NA}
bio ncov --fasta | head -5
```

## Extract a partial sequence and change the sequence id

```{bash, comment=NA}
bio ncov --fasta --start 100 --end 200 --seqid foo | head -5
```

## Extract the last ten DNA bases

```{bash, comment=NA}
bio ncov --fasta --start -10 --seqid last | head -5
```

## Extract the sequences for features of a certain type

```{bash, comment=NA}
bio ncov --fasta --type CDS | head -5
```

## Extract sequence by gene name

```{bash, comment=NA}
bio ncov:S --fasta --start 100 --end 200 | head -5
```

## Extract sequence by feature accession number

```{bash, comment=NA}
bio ncov:YP_009724390.1 --fasta --start 100 --end 200 | head -5
```

## Translate the sequence

This command translates the DNA sequence to peptides:

```{bash, comment=NA}
bio ncov:S --fasta --start 100 --end 200 --translate | head -5
```

The slice 100 to 200 is applied on the DNA sequence before the translation.

## Extract the protein sequence

This flag extracts the protein sequence embedded in the original GenBank file:

```{bash, comment=NA}
bio ncov:S --fasta --start 100 --end 200 --protein | head -5
```

Note how in this case the slice 100 to 200 is applied on the protein sequence.