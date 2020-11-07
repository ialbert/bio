# bio: convert {#convert}

The `bio` package may be used to 

- convert data from GenBank to other formats: FASTA, GFF
- extract only certain sections of data: features by name, coordinate, range

### View the JSON representation of the data:

```{bash, comment=NA}
bio ncov | head
```

### View the JSON representation for a feature:
  
```{bash, comment=NA}
bio ncov --type CDS --gene S | head 
```

### Convert all features to GFF:

```{bash, comment=NA}
bio ncov --gff | head -5
```

### Convert to GFF only the features with type `CDS`

```{bash, comment=NA}
bio ncov --gff --type CDS | head -5
```

### Convert to GFF only the features tagged with gene `S`

```{bash, comment=NA}
bio ncov --gff --gene S | head -5
```

### Convert to GFF only the features that overlap a interval

```{bash, comment=NA}
bio ncov --gff --start 2000 --end 3000 | head -5
```

### Convert the data to FASTA (the origin of the GenBank)

```{bash, comment=NA}
bio ncov --fasta | head -5
```

### Extract a partial sequence and change the sequence id

```{bash, comment=NA}
bio ncov --fasta --start 100 --end 200 --seqid foo | head -5
```

### Extract the last ten DNA bases

```{bash, comment=NA}
bio ncov --fasta --start -10 --seqid last | head -5
```


### Extract the sequences for features of a certain type

```{bash, comment=NA}
bio ncov --fasta --type CDS | head -5
```


### Gene name may be encoded in the name

```{bash, comment=NA}
bio ncov:S --fasta --start 100 --end 200 | head -5
```
