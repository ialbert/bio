# Run information {#bio-runinfo}

`bio` has support to automatically query your data for additional information at sra.

While Not all GenBank records are properly cross referenced, but for those that are `bio` can get you the SRA inforamation right away. Here is how it works. Get a strain of Ebola sequenced in 2014

```{bash, comment=NA}
bio fetch KM233118 --rename ebola 
```

First check to see if the record is being properly cross referenced:

```{bash, comment=NA}
bio runinfo ebola
```

We can see that the data has both a BioProject and a BioSample associated with it. It means we may obtained a more the detailed information on the sequencing data that produced the information:

```{bash, comment=NA}
bio runinfo ebola --sample | head -38
```

We can also obtain the full run information for the entire project (we are limiting the results to make the query speedier):

```{bash, comment=NA}
bio runinfo ebola --project --limit 10 | head -10
```

You can also produce the output in a tab delimited format:

```{bash, comment=NA}
bio runinfo ebola --project --table --limit 10 | cut -f 1,5,8,12,15,19,29 | head
```