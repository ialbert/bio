# Welcome to `bio`

> The software is currently under development. It is operational but not fully vetted.

`bio` - command-line utilities to make bioinformatics explorations more enjoyable. `bio` streamlines the tedious bioinformatics and lets users quickly answer questions such as:
 
- *How do I automate the sequence for a viral genome?*
- *How do I obtain the biological annotation of data?*
- *How do I get the coding sequence for a specific gene?* 
- *What are the differences between two sequences?*
- *What is the lineage of SARS-COV-2?*
- *What are minisatellites and  microsatellites?*
- *What is a "tolerance induction to tumor cell"?* 

`bio` combines and represents data from different sources: GenBank, SRA, Gene Ontology, Sequence Ontology, NCBI Taxonomy through a unified interface. Having access to all the utility described above makes the `bio` package well suited for exploratory analysis of genomes. 

The software was written to teach bioinformatics and is the companion software to the [Biostar Handbook][handbook]
 
[biopython]: https://biopython.org/
[emboss]: http://emboss.sourceforge.net/
[simplesam]: https://github.com/mdshw5/simplesam 
[handbook]: https://www.biostarhandbook.com/
 
## Quick links

* Source code: https://github.com/ialbert/bio
* Documentation: https://www.bioinfo.help

[usage]: https://github.com/ialbert/bio/blob/master/test/bio_examples.sh

## Why does this software exist?

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, reading documentation, and numerous other preparations that slow down your progress. 

Time and again, I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package is meant to solve that tedium. 

## Diving right in

Here is how to align the sequences of SARS-COV-2 (`NC_045512`) versus the same region of a bat coronavirus (`MN996532`):

    # Obtain the data.
    bio NC_045512 MN996532 --fetch 
    
Align the sequences (showing 60bp for brevity).

```{bash, comment=NA}
bio NC_045512 MN996532 --align --end 60  
```

## A more realistic workflow

Suppose you wanted to identify the mutations between the `S` protein of the bat coronavirus deposited as `MN996532` and the `S` protein of the ancestral SARS-COV-2 virus designated by the NCBI via accession number `NC_045512`. 

If you are a trained bioinformatician, think about all the steps you would need to perform to accomplish this task, then think about the effort it would take you to teach someone else how to do the same. 

With the `bio` package, the process takes simple, concise steps.

## Download and rename

First, we download and rename the data keep our sanity:

    bio NC_045512 --fetch --rename ncov
    bio MN996532  --fetch --rename ratg13

From now on, `bio` can operate on  `NC_045512` using the name `ncov` and on `MN996532` using the name `ratg13` no matter where you are on your computer! 

## Convert to different formats

`bio` stores data in an internal storage system that it can find from any location. There is no clutter of files or paths to remember. For example, in any directory, you now can type:

```{bash, comment=NA}
    bio ncov --fasta --end 100 | head -2
```
    
and it will show you the FASTA representation of  the genome     

You could also convert the data stored under `ncov` name to other formats. Let's convert features with type `CDS` to `GFF`:

```{bash, comment=NA}
bio ncov --gff --type CDS  | head -5
```

## Align nucleotides or peptides

Now, back to our problem of aligning proteins. Let's align the first 90 base pairs of DNA sequences for the `S` protein for each organism, `bio` even gives you a shortcut; instead of typing `--gene S --type CDS` you can write it as `ncov:S` :

```{bash, comment=NA}
bio ncov:gene:S ratg13:S --end 60 --align
```
    
We can visualize the translation of the DNA into aminoacids with one letter (`-1`) or three-letter codes (`-3`):  
   
```{bash, comment=NA}
bio ncov:gene:S ratg13:gene:S --end 60 --align -1
```
    
If, instead, we wanted to align the 60bp DNA subsequences for `S` protein after their translation into proteins, we could do it like so:

```{bash, comment=NA}
bio ncov:gene:S ratg13:gene:S --translate --end 60 --align
```
    
We can note right away that all differences in the first 60bp of DNA are synonymous substitutions, the protein translations are the same.


## Look up the taxonomy

`bio` understands taxonomies. Finding the lineage of the organism is as simple as:

```{bash, comment=NA}
bio ncov --taxon --lineage
```

## See the bioproject

`bio` knows about bio projects and sequencing data.
As it turns out the data for `ncov`  data is not adequately cross-referenced at NCBI ... thus, we can't quite get the SRR run numbers automatically, even at NCBI.

Let's pick another data that has better cross-references, perhaps a virus from the 2014 Ebola outbreak:

    bio  KM233118 --fetch --rename ebola14

and now print:

```{bash, comment=NA}
bio ebola14 --sra 
```
   
if we wanted the SRR run numbers, we could run:

    bio ebola14 --sra --sample
 
to get:

     
    [
        {
            "Run": "SRR1553609",
            "ReleaseDate": "2014-08-19 11:41:53",
            "LoadDate": "2014-08-19 11:18:49",
            "spots": "464802",
            "bases": "93890004",
            "spots_with_mates": "464802",
            "avgLength": "202",
            "size_MB": "51",
            "download_path": "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1553609/SRR1553609.1",
            "Experiment": "SRX674271",
            "LibraryName": "NM042.3.FCH9",
            "LibraryStrategy": "RNA-Seq",
            "LibrarySelection": "cDNA",
            "LibrarySource": "TRANSCRIPTOMIC",
            "LibraryLayout": "PAIRED",
    ...
                

## `bio` is a data model

Beyond the functionality that we show, `bio` is also an exploration into modeling biological data. The current standards and practices are woefully antiquated and painfully inadequate. Default formats such as GenBank or EMBL are inefficient and tedious to program with. 

In contrast `bio` represents data in simple, efficient, compressed in JSON format. 

```{bash, comment=NA}
bio ncov | head -20
```

The data layout allows `bio` to read in the entire human chromosome 1, with its 253 million characters and 328 thousand genomic features, in just three(!) seconds. In another 3 seconds, `bio`  can convert that information FASTA or GFF; it can filter it by type, translate the sequence, extract proteins, slice by coordinate, etc.:

    time bio chr1 --fasta | wc -c
    253105766

    real    0m6.238s
    user    0m4.156s
    sys     0m2.172s

For shorter genomes, bacterial or viral, the conversion times are under a fraction of a second.  

Thanks to the representation, it is trivially easy to extend `bio`. The data is already structured in an efficient layout that needs no additional parsing to load. 