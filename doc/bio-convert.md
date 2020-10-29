# bio: conversion

The `bio` package may be used to 

- convert data from GenBank to other formats: FASTA, GFF
- extract only certain sections of data: features by name, coordinate, range

```{r, child='bio-tips.md'}
```

### View the JSON representation of the data:

```{bash, comment=NA}
bio NC_045512 | head
```

### View the JSON representation for a feature:
  
```{bash, comment=NA}
bio NC_045512 --type CDS --gene S | head 
```

### Convert all features to GFF:

```{bash, comment=NA}
bio NC_045512 --gff | head -5
```

### Convert to GFF only the features with type `CDS`

```{bash, comment=NA}
bio NC_045512 --gff --type CDS | head -5
```

### Convert to GFF only the features tagged with gene `S`

```{bash, comment=NA}
bio NC_045512 --gff --gene S | head -5
```

### Convert to GFF only the features that overlap a interval

```{bash, comment=NA}
bio NC_045512 --gff --start 2000 --end 3000 | head -5
```

### Convert the data to FASTA (the origin of the GenBank)

```{bash, comment=NA}
bio NC_045512 --fasta | head -5
```

### Extract a partial sequence and change the sequence id

```{bash, comment=NA}
bio NC_045512 --fasta --start 100 --end 200 --seqid foo | head -5
```

### Extract the last ten DNA bases

```{bash, comment=NA}
bio NC_045512 --fasta --start -10 --seqid last | head -5
```


### Extract the sequences for features of a certain type

```{bash, comment=NA}
bio NC_045512 --fasta --type CDS | head -5
```
