---
title: "The bio package"
documentclass: book
fontsize: 12pt
numbering:  "false"
---
# Welcome to `bio`

`bio` - command-line utilities to make bioinformatics explorations more enjoyable. 

`bio` streamlines the tedious bioinformatics and lets users quickly answer questions such as:
 
- *How do I access a sequence for a viral genome?*
- *How do I obtain the biological annotation of data?*
- *How do I get the coding sequence for a specific gene?* 
- *What are the differences between two sequences?*
- *What is the lineage of SARS-COV-2?*
- *What are minisatellites and  microsatellites?*

`bio` combines and represents data from different sources: [GenBank][genbank], [Gene Ontology][go], [Sequence Ontology][so], 
[NCBI Taxonomy][taxonomy] and [Short Read Archive][sra] through a unified interface. Having access to all the utility described above makes the `bio` package well suited for exploratory analysis of genomes. 

The software was written to teach bioinformatics and is the companion software to the [Biostar Handbook][handbook]
 
[biopython]: https://biopython.org/
[emboss]: http://emboss.sourceforge.net/
[simplesam]: https://github.com/mdshw5/simplesam 
[handbook]: https://www.biostarhandbook.com/
[genbank]: https://www.ncbi.nlm.nih.gov/genbank/
[sra]: https://www.ncbi.nlm.nih.gov/sra
[taxonomy]: https://www.ncbi.nlm.nih.gov/taxonomy
[so]: http://www.sequenceontology.org/
[go]: http://geneontology.org/

## Quick links

* Source code: https://github.com/ialbert/bio
* Documentation: https://www.bioinfo.help

[usage]: https://github.com/ialbert/bio/blob/master/test/bio_examples.sh

## Why does this software exist?

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, reading documentation, and numerous other preparations that slow down your progress. 

Time and again, I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package is meant to solve that tedium. 

## Diving right in

Here is how to align the sequences of SARS-COV-2 (`NC_045512`) versus the same region of a bat coronavirus (`MN996532`). First get the data:

```{bash, eval=FALSE, code=readLines("code/intro-01.sh")}
```

```{r, eval=FALSE, code=readLines("code/intro-01.sh.txt")}
```


Now align the sequences (showing 60bp for brevity).

```{bash, eval=FALSE, code=readLines("code/intro-02.sh")}
```

```{r, eval=FALSE, code=readLines("code/intro-02.sh.txt")}
```


## A more realistic workflow

Suppose you wanted to identify the mutations between the `S` protein of the bat coronavirus deposited as `MN996532` and the `S` protein of the ancestral SARS-COV-2 virus designated by the NCBI via accession number `NC_045512`. 

If you are a trained bioinformatician, think about all the steps you would need to perform to accomplish this task, then think about the effort it would take you to teach someone else how to do the same. 

With the `bio` package, the process takes simple, concise steps.

## Download and rename

First, we download and rename the data keep our sanity:


    bio fetch NC_045512 --rename ncov
    bio fetch MN996532  --rename ratg13


From now on, `bio` can operate on  `NC_045512` using the name `ncov` and on `MN996532` using the name `ratg13` no matter where you are on your computer! 

## Convert to different formats

`bio` stores data in an internal storage system that it can find from any location. There is no clutter of files or paths to remember. For example, in any directory, you now can type:


    bio convert ncov --fasta --end 100 | head -2

    
and it will show you the FASTA representation of  the genome     

You could also convert the data stored under `ncov` name to other formats. Let's convert features with type `CDS` to `GFF`:

    bio convert ncov --gff --type CDS  | head -5

## Align nucleotides or peptides

Now, back to our problem of aligning proteins. Let's align the first 90 base pairs of DNA sequences for the `S` protein for each organism, `bio` even gives you a shortcut; instead of typing `--gene S --type CDS` you can write it as `ncov:S` :

    bio align ncov:S ratg13:S --end 60

We can visualize the translation of the DNA into aminoacids with one letter (`-1`) or three-letter codes (`-3`):  
   
    bio align ncov:S ratg13:S --end 60 -1

If, instead, we wanted to align the 60bp DNA subsequences for `S` protein after their translation into proteins, we could do it like so:


    bio align ncov:S ratg13:S --translate --end 60

    
We can note right away that all differences in the first 60bp of DNA are synonymous substitutions, the protein translations are the same.


## Look up the taxonomy

`bio` understands taxonomies. Finding the lineage of the organism is as simple as:


    bio taxon ncov --lineage


## `bio` is a data model

Beyond the functionality that we show, `bio` is also an exploration into modeling biological data. The current standards and practices are woefully antiquated and painfully inadequate. Default formats such as GenBank or EMBL are inefficient and tedious to program with. 

In contrast `bio` represents data in simple, efficient, compressed in JSON format. 

    bio convert ncov --json | head -20

The data layout allows `bio` to read in the entire human chromosome 1, with its 253 million characters and 328 thousand genomic features, in just three(!) seconds. In another 3 seconds, `bio`  can convert that information FASTA or GFF; it can filter it by type, translate the sequence, extract proteins, slice by coordinate, etc.:

    time bio convert chr1 --fasta | wc -c
    253105766

    real    0m6.238s
    user    0m4.156s
    sys     0m2.172s

For shorter genomes, bacterial or viral, the conversion times are under a fraction of a second.  

Thanks to the representation, it is trivially easy to extend `bio`. The data is already structured in an efficient layout that needs no additional parsing to load. 
