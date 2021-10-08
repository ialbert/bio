# bio align: perform alignments {#bio-align}

> **Note:** alignments in `bio` are primarily designed for exploratory use, for aligning relatively short (up to ~30Kb long sequences), visually investigating the alignments, interacting with the sequences before and after alignment.
>
> Use software that relies on heuristics when investigating large datasets. Specialized software will operate (many) orders of magnitude faster. Depending on your needs you may want to use: `blast`, `blat`, `mummer`, `minimap2`, `lastz`, `lastal`, `exonerate`, `vsearch`, `diamond`.

Install `bio` with:

    pip install bio --upgrade

The full documentation for `bio` is maintained at <https://www.bioinfo.help/>.

## Quick start

Input may be given from command line:

    bio align GATTACA GATCA

the first sequence is the target, all following sequences are aligned considered queries. The command prints:

```{r, code=xfun::read_utf8('code/align1.txt'), eval=F}
```

## Alignment as a table

When generating alignment providing output in different formats is important. Here is the alignment formatted as a table:

    bio align GATTACA GATCA --table

prints:

```{r, code=xfun::read_utf8('code/align2.txt'), eval=F}
```

## Alignment as VCF

    bio align GATTACA GATCA --vcf

prints:


```{r, code=xfun::read_utf8('code/align3.txt'), eval=F}
```

## Alignment as differences

The output may be formatted in as differences (mutations) between the reference (second sequence) and the query (first sequence):

    bio align GATTACA GATCA --diff

prints:

```{r, code=xfun::read_utf8('code/align4.txt'), eval=F}
```

You can think of the mutations output can be thought of as a simplified VCF.

## Alignment types

The default alignment is semi-global (global alignment with no end gap penalies). To perform a global alignment write:

    bio align GATTACA GATCA --global

the output is:

```{r, code=xfun::read_utf8('code/align5.txt'), eval=F}
```

Available ptions are : `--global`, `--local`, `--semi-global`

## Align realistic data

Fetches two genomic files:

    bio fetch MN996532 NC_045512 > genomes.gb

Align two genomes, the first is the target, all following sequences are considered queries and are aligned against the target. Here we are looking at what mutations would need to be made to the bat genome to turn it into the coronavirus genome:

    cat genomes.gb | bio fasta --genome | bio align | head

the command aligns two 30KB sequences and takes about 15 seconds on my system, it will print:

    # MN996532.2 (29855) vs NC_045512.2 (29903)
    # pident=96.0% len=29908 ident=28725 mis=1125 del=5 ins=53
    # semiglobal: score=139005.0 gap open=11 extend=1  matrix=NUC.4.4
    
    ATTAAAGGTTTATACCTTTCCAGGTAACAAACCAACGAACTCTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAA
    ||||||||||||||||||.|||||||||||||||||.||||.|||||||||||||||||||||||||||||||||||||||
    ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAA

The `bio align` method takes the first file that it sees as target and aligns all other sequences to it as queries.

## Align coding sequences

Align the DNA corresponding to protein `S`

    cat genomes.gb | bio fasta --gene S | bio align | head

prints:

```{r, code=xfun::read_utf8('code/align6.txt'), eval=F}
```

## Align protein sequences

Align the protein sequences that prints:

    cat genomes.gb | bio fasta --gene S --protein |  head

that prints:

```{r, code=xfun::read_utf8('code/align7.txt'), eval=F}
```

## Alignment showing mutations

    cat genomes.gb | bio fasta --gene S --protein | bio align --diff | tail -5

prints the variations:

```{r, code=xfun::read_utf8('code/align9.txt'), eval=F}
```

## Alignment with tabular output

You can produce a column based table output

    cat genomes.gb | bio fasta --gene S --protein | bio align --table

prints:

```{r, code=xfun::read_utf8('code/align8.txt'), eval=F}
```


## Different scoring matrices

    cat genomes.gb | bio fasta --gene S --translate | bio align --matrix PAM30 | head

prints:

```{r, code=xfun::read_utf8('code/align10.txt'), eval=F}
```

### Scoring matrices

The scoring matrix may be a builtin name or a file with a scoring matrix format. See the scoring with:

    bio align --matrix PAM30

prints:

    #
    # This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
    #
    # PAM 30 substitution matrix, scale = ln(2)/2 = 0.346574
    #
    # Expected score = -5.06, Entropy = 2.57 bits
    #
    # Lowest score = -17, Highest score = 13
    #
        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    A   6  -7  -4  -3  -6  -4  -2  -2  -7  -5  -6  -7  -5  -8  -2   0  -1 -13  -8  -2  -3  -3  -3 -17
    R  -7   8  -6 -10  -8  -2  -9  -9  -2  -5  -8   0  -4  -9  -4  -3  -6  -2 -10  -8  -7  -4  -6 -17
    N  -4  -6   8   2 -11  -3  -2  -3   0  -5  -7  -1  -9  -9  -6   0  -2  -8  -4  -8   6  -3  -3 -17
    D  -3 -10   2   8 -14  -2   2  -3  -4  -7 -12  -4 -11 -15  -8  -4  -5 -15 -11  -8   6   1  -5 -17
    ...

## Exercises

### Align genomic DNA to CDNA

    bio fetch ENST00000288602 | head > genomic.fa

    bio fetch ENST00000288602 --type cdna | head > cdna.fa

Is this a good alignment?

    bio align genomic.fa cdna.fa

Try local alignment

    bio align genomic.fa cdna.fa --local
