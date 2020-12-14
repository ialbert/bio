# The JSON format {#bio-json}

`bio` obtains data from NCBI and transforms it into an internal, simpler format. One would only need to process this format to provide functionality that is not yet offered in `bio` 

### Get a dataset

Get SARS-COV-2 data and rename it to `ncov`:

```{bash, comment=NA}
bio NC_045512 --fetch --rename ncov
```

### The GenBank data

Explore the contents of the file downloaded from NCBI

```{bash, comment=NA}
bio ncov --genbank | head -20
```

### JSON data representation

See the transformed GenBank file as the JSON representation:

```{bash, comment=NA}
bio ncov --json | head -36
```

## References

The following references may be consulted to understand how data should be represented in GenBank and GFF formats:

INSDC feature descriptions:

* http://www.insdc.org/files/feature_table.html#2

NCBI GenBank format:

* https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

NCBI GFF format:

* https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
