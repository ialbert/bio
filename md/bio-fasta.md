# Convert to FASTA {#bio-fasta}


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

## Gene name may be encoded in the name

```{bash, comment=NA}
bio ncov:S --fasta --start 100 --end 200 | head -5
```
