# Related software {#bio-related}

Since the beginning of bioinformatic times mankind was held down, hindered and obstructed by suboptimal dataformats, clunky interfaces, tedious procedures, bureaucratic protocols. No wonder a thriving cottage industry was born to deal with the challenges, `bio` being one of the proud entrants in the fray.

Here is a list of packages we've evaluated and sometimes borrowed ideas from. Alas none quite works as we needed it, hence like our father, and our father's father and their great grandfather, we have also set out to *finally* fix the problems once and for all. 

Descriptions are direct quotes from the project's homepage.
 
## EMBOSS

EMBOSS is "The European Molecular Biology Open Software Suite". EMBOSS is a free Open Source software analysis package specially developed for the needs of the molecular biology (e.g. EMBnet) user community. The software automatically copes with data in a variety of formats and even allows transparent retrieval of sequence data from the web. Also, as extensive libraries are provided with the package, it is a platform to allow other scientists to develop and release software in true open source spirit. EMBOSS also integrates a range of currently available packages and tools for sequence analysis into a seamless whole. EMBOSS breaks the historical trend towards commercial software packages.

* http://emboss.sourceforge.net/

## Entrez Direct

Entrez Direct (EDirect) provides access to the NCBI's suite of interconnected databases (publication, sequence, structure, gene, variation, expression, etc.) from a Unix terminal window. Search terms are entered as command-line arguments. Individual operations are connected with Unix pipes to allow construction of multi-step queries. Selected records can then be retrieved in a variety of formats.

* https://www.ncbi.nlm.nih.gov/books/NBK179288/

## gffutils

gffutils is a Python package for working with and manipulating the GFF and GTF format files typically used for genomic annotations. Files are loaded into a sqlite3 database, allowing much more complex manipulation of hierarchical features (e.g., genes, transcripts, and exons) than is possible with plain-text methods alone.

* https://github.com/daler/gffutils

## Biocode

This is a collection of bioinformatics scripts many have found useful and code modules which make writing new ones a lot faster.

* https://github.com/jorvis/biocode

## GenomePy

Easily install and use genomes in Python and elsewhere!

The goal is to have a simple and straightforward way to download and use genomic sequences. Currently, genomepy supports UCSC, Ensembl and NCBI.

* https://github.com/vanheeringen-lab/genomepy

## NCBI genome download

Some script to download bacterial and fungal genomes from NCBI after they restructured their FTP a while ago.

Idea shamelessly stolen from Mick Watson's Kraken downloader scripts that can also be found in Mick's GitHub repo. However, Mick's scripts are written in Perl specific to actually building a Kraken database (as advertised).

So this is a set of scripts that focuses on the actual genome downloading.

* https://github.com/kblin/ncbi-genome-download


## Comparisons to EMBOSS

The software with the most similar goals to `bio` is the [emboss suite][emboss], a revolutionary software package developed decades ahead of its time. Unfortunately, perhaps because of being developed so early on, the amazing feats of software engineering within `emboss` are deployed with a nearly incomprehensible documentation that attempts, in vain, to describe an incredibly obtuse command interface. 

We love the concept of `emboss` but even after many years we don't understand how to use it. We constantly have to consult the manual for details. Moreover commands that use `emboss` suites tend to end up as a series of hard to read arcane commands that are surprisingly difficult to comprehend even for experienced scientists. 

Criticism aside, imitation is the greatest form of flattery, `bio` is an homage to `emboss` with the hope that one day, we can replace some functionality from `emboss` with code that brings joy rather than frustrations. 

