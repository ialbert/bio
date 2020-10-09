# bio: view

The `view` command may be used to:

- convert data from GenBank/EMBL to other formats: FASTA, GFF
- extract sections of data: features by name, coordinate, range

## Usage

```{bash, comment=NA}
bio view -h
```

## Examples

The following sections show a command and its output. For clarity we truncate the output to first few relevant lines.

If  add multiple accession numbers the operations will take place sequentially on each.

You may also combine multiple parameters, in that case each condition will be applied.

Coordinates are 1 based (inclusive on both ends) identical to GFF coordinate formats.

### View the GenBank file:

```{bash, comment=NA}
bio view NC_045512 | head
```

### Convert all features to GFF:

```{bash, comment=NA}
bio view NC_045512 --gff | head -3
```

### Convert only features with type `CDS`

```{bash, comment=NA}
bio view NC_045512 --gff --type CDS | head -3
```

### Convert only features that related to gene `S`:

```{bash, comment=NA}
bio view NC_045512 --gff --gene S | head -3
```

### Convert features that overlap a range:

```{bash, comment=NA}
bio view NC_045512 --gff --start 21563 --end 21565 | head -3
```

### Convert the "origin" of the GenBank to FASTA:

```{bash, comment=NA}
bio view NC_045512 --fasta | head -5
```

### Extract a subsequence from the "origin" to FASTA:

```{bash, comment=NA}
bio view NC_045512 --fasta -start 100 -end 200 | head -5
```
