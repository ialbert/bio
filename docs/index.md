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

## Quickstart example

Suppose you wanted to align the sequences of SARS-COV-2 (`NC_045512`) versus the same region of a bat coronavirus (`MN996532`).

### Fetch data

This is how to download the data so that `bio` can operate on it:

```{bash, child='code/index-fetch.txt'}
```

### Align the genomes

Now align the sequences (showing 60bp for brevity).

```{bash, child='code/index-align.txt'}
```

### Convert to FASTA format

Bioinformatics workflows require you to have data in different formats. `bio` can convert data for you.

```{bash, child='code/index-fasta.txt'}
```

### Convert to GFF format:

```{bash, child='code/index-gff.txt'}
```

View the resulting files in IGV



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
