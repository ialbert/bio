## JSON representation {#json}

`bio` obtains data from NCBI and transforms it into an internal, simpler format.

You would only need to process this format to provide functionality that is not yet offered in `bio` 

### Get a dataset

Get SARS-COV-2 data and rename it to `ncov`:

```{bash, comment=NA}
    bio NC_045512 --fetch --rename ncov
```

### The GenBank content

Explore the contents of the file downloaded from NCBI

```{bash, comment=NA}
bio ncov --genbank | head 
```

### JSON data representation

See the transformed GenBank file as the JSON representation:

```{bash, comment=NA}
bio ncov | head
```

### The JSON for a feature
  
Filter the internal representation for a type and a gene name:
  
```{bash, comment=NA}
bio ncov --type CDS --gene S | head 
```

## References

The following references may be consulted to understand how data should be represented in GenBank and GFF formats:

INSDC feature descriptions:

* http://www.insdc.org/files/feature_table.html#2

NCBI GenBank format:

* https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

NCBI GFF format:

* https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
