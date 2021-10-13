# bio table: convert to table {#bio-table}

Converts input to tabular formats. Suppose you get the following data:

    bio fetch NC_045512 MN996532 > genomes.gb

To turn the GenBank into a table run:

    cat genomes.gb | bio table | head

by default it prints the id, size associated with each sequence:

```{r, code=xfun::read_utf8('code/table1.txt'), eval=F}
```

a number of additional fields may be specified, for example:

    cat genomes.gb | bio table --fields id,type,size,gene | head

connects protein ids to genes:

```{r, code=xfun::read_utf8('code/table2.txt'), eval=F}
```

## Additional fields

For a simple enumeration of of the sequence ids and their sizes:

    cat genomes.gb | bio table --fields id,size

or to extract additional metadata:

    cat genomes.gb | bio table --type CDS --fields id,gene,isolate,country,date | head

prints:

```{r, code=xfun::read_utf8('code/table3.txt'), eval=F}
```
