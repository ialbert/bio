# bio format: reformat alignments {#bio-format}

Note: the inputs for this tool must be in the so called "aligned" FASTA. The sequences need to arranged to contain gaps in the places of insertions and deletions. For example:


    >seq1
    THIS-LI-NE-
    >seq2
    --ISALIGNED

Many tools can generate such alignments.

    mamba install mafft

## Usage example

Get FASTA sequences for the coronavirus and the closest related bat coronavirus

    bio fetch MN996532 NC_045512 | bio fasta --genome > genomes.fa

Align the fasta files with `mafft`:

    mafft --auto --quiet --preservecase genomes.fa  > aligned.fa

## Pairwise alignment

The default output generates the diffs:

    cat aligned.fa | bio format | head -12

shows the pairwise alignment:

```{r, code=xfun::read_utf8('code/format1.txt'), eval=F}
```

## Show differences

    cat aligned.fa | bio format --mut | head -5

This output generates an easy-to-read output of the alignments:

```{r, code=xfun::read_utf8('code/format2.txt'), eval=F}
```

Columns are: position, type, query, change, target name.

## Generate VCF output

    cat aligned.fa | bio format --vcf | head

prints:

```{r, code=xfun::read_utf8('code/format3.txt'), eval=F}
```


