---
title: "The bio package"
documentclass: book
fontsize: 12pt
numbering:  "false"
---
# Welcome to `bio`

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, reading documentation, and other preparations that slow down progress. Time and again, I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package was designed  to solve that tedium.

The `bio` software was written to make bioinformatics explorations more enjoyable. The software lets users quickly answer questions such as:
 
- *How do I access a sequence for a viral genome?*
- *How do I obtain the biological annotation of data?*
- *How do I get the coding sequence for a specific gene?*
- *What is the lineage of SARS-COV-2?*
- *What are minisatellites and  microsatellites?*

`bio` combines data from different sources: [GenBank][genbank], [Gene Ontology][go], [Sequence Ontology][so],
[NCBI Taxonomy][taxonomy] and provides an logical, unified interface. The utility described above makes the `bio` package exceedingly well suited for exploratory analysis of genomes.

The software is also used to demonstrate and teach bioinformatics and is the companion software to the [Biostar Handbook][handbook].
 
[biopython]: https://biopython.org/
[emboss]: http://emboss.sourceforge.net/
[simplesam]: https://github.com/mdshw5/simplesam 
[handbook]: https://www.biostarhandbook.com/
[genbank]: https://www.ncbi.nlm.nih.gov/genbank/
[sra]: https://www.ncbi.nlm.nih.gov/sra
[taxonomy]: https://www.ncbi.nlm.nih.gov/taxonomy
[so]: http://www.sequenceontology.org/
[go]: http://geneontology.org/


[usage]: https://github.com/ialbert/bio/blob/master/test/bio_examples.sh

## Quickstart

Suppose you wanted to align the sequences of SARS-COV-2 (`NC_045512`) versus the same region of a bat coronavirus (`MN996532`). Here is how it would work with `bio`.

### 1. Fetch data

Ddownload the data so that `bio` can operate on it. This step needs to be peformed once only:

```{bash, child='code/index-fetch.txt'}
```

### 2. Align the genomes

Now align the sequences (showing just 60bp for brevity).

```{bash, child='code/index-align.txt'}
```

## Data conversion

Bioinformatics workflows often requires you to present data in different formats. `bio` can convert it for you on the fly:

### Convert to FASTA format

```{bash, child='code/index-fasta.txt'}
```

### Convert to GFF format:

```{bash, child='code/index-gff.txt'}
```

View the resulting files in IGV

```{r fig.align='center', echo=FALSE}
knitr::include_graphics('images/igv-index.png', dpi = NA)
```

Among the many useful features, `bio` is also able to generate beautiful data models from GenBank file.

## Data integration

`bio` understands taxonomies, NCBI bioprojects and metadata.

### Navigate the taxonomy

 Finding the lineage of the organism in a GenBank file is as simple as:

```{bash, child='code/index-taxon.txt'}
```

### Getting sample metadata

Get sample metadata for the viral genomes (taxid `2697049`):

```{bash, child='code/index-meta.txt'}
```

## Where to go next

Look at the sidebar for detailed documentation on how `bio` operates.

## Other links

* Source code: https://github.com/ialbert/bio
* Documentation: https://www.bioinfo.help
