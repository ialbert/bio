# bio: introduction

> Beta testing (in progress)

**Making bioinformatics fun again.**

Command-line utilities to make bioinformatics explorations more enjoyable.

Built on top of [BioPython][bioython] and other existing packages; the `bio` software package streamlines bioinformatics tasks such as:
 
- downloading data from NCBI
- converting between data formats 
- extracting information from files (by gene, by coordinate etc)
- aligning sequences

The `bio` package is well suited for exploratory analysis of genomes. 

[biopython]: https://biopython.org/

## Rationale

If you've ever done bioinformatics you know how even seemingly straigthforward tasks require multiple steps, arcane incantation, reading obtuse documentation and several preparations that slow down progress. Time and again I found myself not pursuing an idea because getting to the fun part was too tedios. `bio` is meant to solve that tedium.

For example, suppose you wanted to identify the differences between the `S` protein of the bat coronavirus deposited as `MN996532` and the `S` protein of the ancestral SARS-COV-2 virus designated by the NCBI via accession number `NC_045512`. If you are a trained bioinformatician, think about all the steps you would need to perform to be succesful at this task, the think about the effort it would take to teach someone else how to do it. 
 
Well, with the `bio` package you can just write:

    bio MN996532 --fetch --rename bat
    bio NC_045512--fetch --rename sars2
    
to get the data, and rename it into more manageable lables. Then you can simply align the whole genomes:

    bio bat sars2 --align 

or align just the genomic region that forms the `S` protein:

    bio bat:S sars2:S --align

What did `bio` do in the backround?
 
1. fetches the data from NCBI
1. creates a more efficient local representation of it
1. stores this representation so that next time no internet connection is necessary
1. generate a global DNA alignment. 

But wait there is more. Perhaps you needed local alignments, no problem:

    bio bat:S sars2:S --align --mode local

or that you wanted to align the sequences as proteins:

    bio bat:S sars2:S --align --protein
   
proteins and local alignments:

    bio bat:S sars2:S --align --protein --mode local
   
proteins and local alignments and just a certain region:
   
    bio bat:S sars2:S --align --protein --mode local --start 100 --end 200 
   
look ma, I can even align the last two aminoacids

    bio bat:S sars2:S --align --protein --start -10 
   
There is a lot more the `bio` than just alignments though. The package is an effort to solve the most nagging and annoying limitations that practitioners have that often boil down to do I extract the data I need and I know is embedded in the file.

## Documentation

The documentation is maintained at

    https://bio.github.io

Or in the github repository as markdown files:

    https://github.com/ialbert/bio/tree/master/doc
