# Welcome to bio {#bio-intro}

If you've ever done bioinformatics you know how even seemingly straightforward tasks often require multiple steps, searching the web, reading documentation, clicking around various websites that all together can slow down your progress.

Time and again, I found myself not pursuing ideas because getting to the fun part was too tedious. The `bio` package was designed  to solve that tedium by making bioinformatics explorations more enjoyable. The software lets users quickly answer questions such as:
 
- How do I download data by an accession number? `bio fetch NC_045512 > genome.gb`
- How do I view the sequence annotation of the data? `bio gff genome.gb`
- How do I get the coding sequence for a specific gene? `bio fasta --gene S genome.gb`
- What is the taxid for SARS-COV-2?  `bio taxon sars2`
- What is the lineage of SARS-COV-2? `bio taxon 2697049 --lineage`
- At what date were viral samples collected? `bio data 2697049 | head`
- What is a  microsatellites? `bio explain microsatellite`
- What genes match the word `HBB`? `bio search HBB --species human`

`bio` combines data from different sources: [GenBank][genbank],[Ensembl][ensembl],[Gene Ontology][go],[Sequence Ontology][so],
[NCBI Taxonomy][taxonomy] and provides an unified, consistent interface.

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
[ensembl]: https://www.ensembl.org

[usage]: https://github.com/ialbert/bio/blob/master/test/bio_examples.sh

## Quick start

Install `bio`:

    pip install bio --upgrade

## Code and info

* Source code: https://github.com/ialbert/bio
* Documentation: https://www.bioinfo.help
* Feedback/Comments/Issues: https://github.com/ialbert/bio/issues


## Usage examples

Suppose we found the accession number to data of interest: `NC_045512` representing the Wuhan-Hu-1 isolate of the coronavirus and
the  `MN996532` that stores information on the most similar bat coronavirus. We would like to investigate the differences in the S protein betweeen the two organisms. Here is how you could do it with `bio`:

### Fetch data

    bio fetch NC_045512 MN996532 > genomes.gb

### Convert to FASTA

    bio fasta genomes.gb  | head

prints:

    >NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
    GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
    CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC
    TTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTT

### Align sequences

From the `genomes.gb` obtained above align the protein sequences for the `S` gene

	cat genomes.gb | bio fasta --gene S --protein | bio align | head

prints:

	# PEP: YP_009724390.1 (1,273) vs QHR63300.2 (1,269) score=6541.0
	# Alignment: length=1273 ident=1240/1273(97.4%) mis=29 del=4 ins=0 gap=4
	# Parameters: matrix=BLOSUM62 gap-open=11 gap-extend=1

	MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDN
	|||||||||||||||||||||||||||||||.|||||||||||||||||.|||||||||||||||||||||||||.||||| 81
	MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSSTRGVYYPDKVFRSSVLHLTQDLFLPFFSNVTWFHAIHVSGTNGIKRFDN

### Show alignment as differences

	cat genomes.gb | bio fasta --gene S --protein | bio align --diff | tail

prints the type and variant at each location:

    493	  SNP	YP_009724390.1	Q/Y	    QHR63300.2
    494	  SNP	YP_009724390.1	S/R	    QHR63300.2
    498	  SNP	YP_009724390.1	Q/Y	    QHR63300.2
    501	  SNP	YP_009724390.1	N/D	    QHR63300.2
    505	  SNP	YP_009724390.1	Y/H	    QHR63300.2
    519	  SNP	YP_009724390.1	H/N	    QHR63300.2
    604	  SNP	YP_009724390.1	T/A	    QHR63300.2
    680	  INS	YP_009724390.1	PRRA/-	QHR63300.2
    1121  SNP	YP_009724390.1	N/S	    QHR63300.2
    1224  SNP	YP_009724390.1	V/I	    QHR63300.2

It shows that at position `680` the coronavirus has an four aminoacid insertion `PRRA` the so called furin-cleavage.

### Visualize the genome data

Convert to GFF

    bio gff genomes.gb  | head

prints:

    ##gff-version 3
    NC_045512.2	.	source	0	29903	.	+	.	ID=1;Name=NC_045512;Parent=NC_045512.2
    NC_045512.2	.	five_prime_UTR	0	265	.	+	.	ID=2;Name=five_prime_UTR-1;Parent=five_prime_UTR-1;color=#cc0e74
    NC_045512.2	.	gene	265	21555	.	+	.	ID=3;Name=ORF1ab;Parent=ORF1ab;color=#cb7a77
    NC_045512.2	.	CDS	265	13468	.	+	.	ID=4;Name=YP_009724389.1;Parent=YP_009724389.1
    NC_045512.2	.	CDS	13467	21555	.	+	.	ID=5;Name=YP_009724389.1;Parent=YP_009724389.1

Load and view the resulting files in IGV

```{r fig.align='center', echo=FALSE}
knitr::include_graphics('images/igv-index.png', dpi = NA)
```

Among the many useful features, `bio` is also able to generate informative gene models from a  GenBank file.

