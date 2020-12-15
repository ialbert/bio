# Convert to GFF {#bio-gff}

## Building a nicer gene model

`bio` creates more meaningful  and nicer GFF visualizations:

    # Get chromosome 2L for Drosophila melanogaster (fruit-fly)
    bio NT_033779 --fetch --rename fly 

convert it to gff:

    bio fly --gff > annotations.gff
    
## GFF created with `bio`
   
Here is a region from the GFF file created with the code above as visualized in IGV:

```{r fig.align='center', echo=FALSE}
knitr::include_graphics('images/gff-model-bio.png', dpi = NA)
```

* Exons will have `transcript_id` and `gene_id` attributes set.
* CDS features have `protein_id` and `gene_id` attributes set.

[SO]: http://www.sequenceontology.org/

### Get a dataset

Get SARS-COV-2 data and rename it to `ncov`:

```{bash, comment=NA}
bio NC_045512 --fetch --rename ncov
```

## Convert all features to GFF:

```{bash, comment=NA}
bio ncov --gff | head -5
```

## Convert to GFF only the features with type `CDS`

```{bash, comment=NA}
bio ncov --gff --type transcript,exon,mRNA,CDS | head -5
```

## Convert to GFF only the features tagged with gene `S`

```{bash, comment=NA}
bio ncov --gff --gene S | head -5
```

