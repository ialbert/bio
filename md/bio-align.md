# Sequence alignments {#align}

The alignments in `bio` are primarily designed for exploratory use, for aligning short (a few thousand basepair long sequences), visually investigating the alignments, interacting with the sequences before and after alignment. In such cases the alignments will be generated within less than a second. The implementations are mathematically optimal but the libraries that we rely on do not scale well to longer sequences.

Use a specially designed software, that relies on heuristics, to perform large scale studies. These tools will operate orders of magnitude faster. Depending on your needs `blast`, `blat`, `mummer`, `minimap2`, `lastz`, `exonerate` will be far better suited for genome wide analyses.
 
## DNA alignment

```{bash, comment=NA}
# Align the extracted protein.
bio ncov:S ratg13:S --end 90 --align
```

## Align the translated regions

```{bash, comment=NA}
bio ncov:S ratg13:S --end 90 --translate --align 
```

## Default alignment is semiglobal

A global alignment where end gaps are have no penality.

```{bash, comment=NA}
bio THISLINE ISALIGNED  -i --align
```

## Local alignment

```{bash, comment=NA}
bio THISLINE ISALIGNED -i --align --local
```

## Global alignment

```{bash, comment=NA}
bio THISLINE ISALIGNED -i --align --global
```

## Semiglobal alignment

```{bash, comment=NA}
bio  THISLINE ISALIGNED -i --align --semiglobal
```
