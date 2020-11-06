# bio: introduction

**Making bioinformatics fun again.**

Command-line utilities to make bioinformatics explorations more enjoyable.

Built on top of [BioPython][biopython] and other existing packages; the `bio` software package streamlines bioinformatics tasks such as:
 
- downloading data from NCBI
- converting between data formats 
- extracting information from files (by gene, by coordinate etc)
- aligning sequences
- and ... many more

Having acces to all the utility above makes the `bio` package well suited for exploratory analysis of genomes. 

The software with the most similar goals is the [emboss suite][emboss].

[biopython]: https://biopython.org/
[emboss]: http://emboss.sourceforge.net/

## Rationale

If you've ever done bioinformatics you know how even seemingly straigthforward tasks require multiple steps, arcane incantations, reading documentation and other preparations that slow down your progress. 

Time and again I found myself not pursuing an idea because getting to the fun part was too tedious. The `bio` package is meant to solve that tedium. 

For example, suppose you wanted to identify the differences between the `S` protein of the bat coronavirus deposited as `MN996532` and the `S` protein of the ancestral SARS-COV-2 virus designated by the NCBI via accession number `NC_045512`. If you are a trained bioinformatician, think about all the steps you would need to perform to accomplish this task, the think about the effort it would take to teach someone else how to do it. Right?
 
Well, with the `bio` package you can just write:

    bio MN996532 --fetch --rename ratg13
    bio NC_045512 --fetch --rename ncov
    
to get the data, and rename it into more manageable labels. Since you are interested in the `S` protein alone you can write:


    bio align ncov:S ratg13:S | head -20

to see:
    
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
    
    QHR63300.2   ATACACCAACTCATCCACCCGTGGTGTCTATTACCCTGACAAAGTTTTCAGATCTTCAGTTTTACATTTAACTCAGGATT
                 ||||||.||.||.|.|||.||||||||.||||||||||||||||||||||||||.|||||||||||||.|||||||||.|
    YP_009724390 ATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGAC

or maybe you wanted to align the `S` protein in protein space like so:

    bio align ncov:S ratg13:S --protein | head -20
    
that prints:

    ### 1: YP_009724390 vs QHR63300.2 ###
    
    Length:	1273 (local) 
    Query:	1273 [1, 1273]
    Target:	1269 [1, 1269]
    Score:	6541
    Ident:	1240/1273 (97.4%)
    Simil:	1252/1273 (98.4%)
    Gaps:	4/1273 (0.3%)
    Matrix:	blosum62(-11, -1) 
    
    QHR63300.2   MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSSTRGVYYPDKVFRSSVLHLTQDLFLPFFSNVTWFHAIHVSGTNGIKRFD
                 |||||||||||||||||||||||||||||||.|||||||||||||||||.|||||||||||||||||||||||||.||||
    YP_009724390 MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFD
    
    QHR63300.2   NPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY
                 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    YP_009724390 NPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY


We can immediately see that many of the mutations are synomous, the similarity went from `92%` to `97%`.

What could `bio` do for us?
 
1. fetch the data from NCBI
1. create a more efficient local representation the data
1. store this representation so that next time no internet connection is necessary
1. generate alignments 

But wait there is more. Perhaps you wanted the FASTA sequence of the genome

    bio ncov --fasta | head -5
    
prints:

    >ncov Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
    GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
    CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC
    TTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTT

But wait there is even more, a lot more. How about translating the, reverse of the last 10 nucleotides of every feature labeled as CDS. Why not:


    bio ncov --fasta --type CDS --start -10 --translate | head -5
    
ah yes, just what I needed:    
        
    >YP_009724389.1 [-9:21291], translated DNA
    *QL
    
    >YP_009725295.1 [-9:13218], translated DNA
    CGV

Or what about GFF regions of type `gene` that overlap with a region 1000 to 2000

    bio ncov --gff --start 1000 --end 2000 | head

it prints:

    ##gff-version 3
    ncov	.	gene	266	21555	.	+	1	Name=ORF1ab;type=gene;gene=ORF1ab;db_xref=GeneID:43740578
    ncov	.	CDS	266	21555	.	+	1	Name=YP_009724389.1;type=CDS;gene=ORF1ab;protein_id=YP_009724389.1;product=ORF1ab polyprotein;db_xref=GeneID:43740578
    ncov	.	mature_protein_region	806	2719	.	+	1	Name=YP_009725298.1;type=mat_peptide;gene=ORF1ab;protein_id=YP_009725298.1;product=nsp2
    ncov	.	CDS	266	13483	.	+	1	Name=YP_009725295.1;type=CDS;gene=ORF1ab;protein_id=YP_009725295.1;product=ORF1a polyprotein;db_xref=GeneID:43740578
    ncov	.	mature_protein_region	806	2719	.	+	1	Name=YP_009742609.1;type=mat_peptide;gene=ORF1ab;protein_id=YP_009742609.1;product=nsp2

## Documentation

The documentation is maintained at

    https://bio.github.io

Or in the github repository as markdown files:

    https://github.com/ialbert/bio/tree/master/doc

## Comparisons to EMBOSS

The software with the most similar goals to `bio` is the [emboss suite][emboss] a package developed way ahead of its time, perhaps the main reason why its amazing feats of software engineering are packaged with incomprehensible documentation and incredibly obtuse user interfaces. 

We love the concept of `emboss` but even after many years we don't fully understand it intricacies, We constantly have to consult the manual for details. Commands that use `emboss` suites always end up as a series of hard to read jumbles of commands that are surprisingly difficult to comprehend even for experienced scientists.

