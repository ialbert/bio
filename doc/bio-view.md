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

You may add multiple accession numbers and the operations will take place on each sequentially.

You may also combine multiple parameters, each condition will be applied.

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
bio view NC_045512 --gff  --type CDS | head -3
```

### Convert only features with name `S`:

```{bash, comment=NA}
bio view NC_045512 --gff  --name S | head -3
```

### Convert features that overlap a range:

```{bash, comment=NA}
bio view NC_045512 --gff  --start 21563 --end 21565 | head -3
```

### Convert the "origin" sequence to FASTA:

```{bash, comment=NA}
bio view NC_045512 --fasta | head -5
```

### Extract DNA sequences for all features



