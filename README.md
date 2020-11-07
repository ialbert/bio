# bio: introduction

**Making bioinformatics fun again.**

Command-line utilities to make bioinformatics explorations more enjoyable.

Built on top of [BioPython][biopython] and other existing packages; the `bio` software package streamlines bioinformatics tasks such as:
 
- downloading data from NCBI
- converting between data formats 
- extracting information from files (by gene, by coordinate etc)
- aligning sequences
- and ... many more

Having acces to all the utility described above makes the `bio` package well suited for exploratory analysis of genomes. 

[biopython]: https://biopython.org/
[emboss]: http://emboss.sourceforge.net/

## Rationale

If you've ever done bioinformatics you know how even seemingly straightforward tasks require multiple steps, arcane incantations, reading documentation and numerous other preparations that slow down your progress. 

Time and again I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package is meant to solve that tedium. 

For example, suppose you wanted to identify the differences between the `S` protein of the bat coronavirus deposited as `MN996532` and the `S` protein of the ancestral SARS-COV-2 virus designated by the NCBI via accession number `NC_045512`. If you are a trained bioinformatician, think about all the steps you would need to perform to accomplish this task, the think about the effort it would take to teach someone else how to do it.
 
Well, with the `bio` package it would work like so. First we get and rename the data to have more manageable labels:

    bio MN996532 --fetch --rename ratg13
    bio NC_045512 --fetch --rename ncov
    
Let's align the first 80 basepairs of DNA sequences for the `S` protein as annotated in each organism:

    bio align ncov:S ratg13:S --end 80

the above command produces:
    
    ### 1: YP_009724390 vs QHR63300.2 ###
    
    Length:	3827 (local) 
    Query:	3822 [1, 3822]
    Target:	3810 [1, 3810]
    Score:	16694
    Ident:	3554/3827 (92.9%)
    Simil:	3554/3827 (92.9%)
    Gaps:	22/3827 (0.6%)
    Matrix:	nuc44(-11, -1) 
    
    QHR63300.2   ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTTTCTAGTCAGTGTGTTAATCTAACAACTAGAACTCAGTTACCTCCTGC
                 ||||||||||||||||||||||||||||||||.||||||||||||||||||||.|||||.||||||||.|||||.|||||
    YP_009724390 ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGC
    
    
If instead we wanted to align the 80bp DNA sequences for `S` protein after their translation into proteins we could do it like so:

    bio align ncov:S ratg13:S --translate --end 80
    
Now the output is:

    ### 1: YP_009724390 vs QHR63300.2 ###
    
    Length: 26 (local)
    Query:  26 [1, 26]
    Target: 26 [1, 26]
    Score:  131
    Ident:  26/26 (100.0%)
    Simil:  26/26 (100.0%)
    Gaps:   0/26 (0.0%)
    Matrix: blosum62(-11, -1)
    
    QHR63300.2   MFVFLVLLPLVSSQCVNLTTRTQLPP
                 ||||||||||||||||||||||||||
    YP_009724390 MFVFLVLLPLVSSQCVNLTTRTQLPP

We can note right away that all differences in DNA are synonymous substitutions, both pieces of DNA code for the same proteins.

What did `bio` do for us?
 
1. fetched the data from NCBI
1. created a more efficient local representation the data
1. stored this representation so that next time no internet connection is necessary
1. generated alignments 

But wait there is more. Perhaps we wanted a subsection of  FASTA sequence of the genome

    bio ncov --fasta --end 160
    
prints:
    
    >ncov Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
    GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
    CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGG

But wait, there is even more, a lot more. How about translating the reverse of the last 10 nucleotides of every feature labeled as `CDS`. `bio` can do that like so:

    bio ncov --fasta --type CDS --start -10 --translate
    
ah yes, just what I needed:    
        
    >YP_009724389.1 [-9:21291], translated DNA
    *QL
    
    >YP_009725295.1 [-9:13218], translated DNA
    CGV
    ...
    
Or what about GFF regions of type `gene` that overlap with a region 1000 to 2000

    bio ncov --gff --start 1000 --end 2000 | head

it prints:

    ##gff-version 3
    ncov	.	gene	266	21555	.	+	1	Name=ORF1ab;type=gene;gene=ORF1ab;db_xref=GeneID:43740578
    ncov	.	CDS	266	21555	.	+	1	Name=YP_009724389.1;type=CDS;gene=ORF1ab;protein_id=YP_009724389.1;product=ORF1ab polyprotein;db_xref=GeneID:43740578
    ncov	.	mature_protein_region	806	2719	.	+	1	Name=YP_009725298.1;type=mat_peptide;gene=ORF1ab;protein_id=YP_009725298.1;product=nsp2
    ncov	.	CDS	266	13483	.	+	1	Name=YP_009725295.1;type=CDS;gene=ORF1ab;protein_id=YP_009725295.1;product=ORF1a polyprotein;db_xref=GeneID:43740578
    ncov	.	mature_protein_region	806	2719	.	+	1	Name=YP_009742609.1;type=mat_peptide;gene=ORF1ab;protein_id=YP_009742609.1;product=nsp2

And so on. `bio` has a wealth of utility that makes bioinformatics more accessible.

## Documentation

The documentation is maintained at

    https://bio.github.io

Or in the github repository as markdown files:

    https://github.com/ialbert/bio/tree/master/doc

## Comparisons to EMBOSS

The software with the most similar goals to `bio` is the [emboss suite][emboss] a revolutionary software package developed decades ahead of its time. Alas perhaps the seminal nature of `emboss` is also the reason why its amazing feats of software engineering are packaged with nearly incomprehensible documentation and uncommonly obtuse user interfaces. 

We love the concept of `emboss` but even after many years we don't understand how to use it. We constantly have to consult the manual for details. Moreover commands that use `emboss` suites tend to end up as a series of hard to read jumbles of commands that are surprisingly difficult to comprehend even for experienced scientists.

`bio` is an homage to `emboss` with the hope that one day we can replace all the functionality from `emboss` in a form that brings joy rather than frustrations.

