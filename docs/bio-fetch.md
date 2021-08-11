# `fetch`: download GenBank data {#bio-fetch}

`bio` can fetch data in GeneBank format based on accession numbers

Multiple ways to fetch data.

You may separate accession numbers with commas

    bio fetch NC_045512,MN996532 > genomes.gb

You may separate accession numbers with spaces

    bio fetch NC_045512 MN996532 > genomes.gb

You may pipe accession numbers into the tool

    echo NC_045512 | bio fetch > genomes.gb
