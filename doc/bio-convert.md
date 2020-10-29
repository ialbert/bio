# bio: convert

The `bio` package may be used to 

- convert data from GenBank to other formats: FASTA, GFF
- operate on certai sections of data: features by name, coordinate, range

## Usage

```
bio -h
```

## Examples

The following sections show a command and its output. For clarity we truncate the output to first few relevant lines.

If  add multiple accession numbers the operations will take place sequentially on each.

You may also combine multiple parameters, in that case each condition will be applied.

Coordinates are 1 based (inclusive on both ends) identical to GFF coordinate formats.

### View the JSON representation of the data:

```{bash, comment=NA}
bio NC_045512 | head
```

### Convert all features to GFF:

```{bash, comment=NA}
bio NC_045512 --gff | head -3
```

### Convert to GFF the features with type `CDS`

```{bash, comment=NA}
bio NC_045512 --gff --type CDS | head -3
```

### Convert to GFF only the features tagged with gene `S`:

```{bash, comment=NA}
bio NC_045512 --gff --gene S | head -3
```

### Convert to GFF features that overlap a range:

```{bash, comment=NA}
bio NC_045512 --gff --start 21563 --end 21565 | head -3
```

### Convert to FASTA the origin sequence in GenBank:

```{bash, comment=NA}
bio NC_045512 --fasta | head -5
```

### Extract a subsequence from the "origin":

```{bash, comment=NA}
bio NC_045512 --fasta --start 100 --end 200 | head -5
```

### Extract the sequences for features of a certain type:

```{bash, comment=NA}
bio NC_045512 --fasta --type CDS | head -5
```
