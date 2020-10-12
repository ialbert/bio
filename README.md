# bio: introduction

> Note: the package is not yet released, some functionality is not yet operational.

**Making learning bioinformatics fun again.**

Command line utilities to make bioinformatics education more streamlined.

Build on top of [BioPython][bioython] the `bio` software package attempts to streamline several bioinformatics tasks,
to allow learners to focus on the concepts that matter.

The primary use of the package is for education and for exploratory analysis of viral and bacterial genomes.

[biopython]: https://biopython.org/

## Rationale

Typical bioinformatics solutions end up  being unnecessarily complicated. Seemingly simple tasks require lengthy command
preparations that slow down progress.

For example suppose you wanted to find the alignment and the differences between protein `S` of the bat corona virus deposited as 
`MN996532` and protein `S` of the ancestral SARS-COV-2 virus (accession numbers as designated by NCBI accession number `NC_045512`). 
If you are a trained bioinformatician think about all the steps you would need to undertake to perform that task.
 
With `bio` you can write:

    bio align MN996532:S NC_045512:S

And that's it. It will:
 
1. automatically fetch the data
1. store the compressed data in a cache so next time it won't need to connect to the internet
1. produce a global DNA alignment. 

But wait there is more. Perhaps you needed local alignments, no problem:

    bio align MN996532:S NC_045512:S --local

or align the translated sequences as proteins:

    bio align MN996532:S NC_045512:S --translation

perhaps you wanted to align only a range of the sequences of protein `S`:

    bio align MN996532:S NC_045512:S --range 100-200

## Documentation

The documentation is maintained at

    https://bio.github.io

Or in the github repository as markdown files:

    https://github.com/ialbert/bio/tree/master/doc

