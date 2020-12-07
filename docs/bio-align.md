# Sequence alignments {#bio-align}

The alignments in `bio` are primarily designed for exploratory use, for aligning short (a few thousand basepair long sequences), visually investigating the alignments, interacting with the sequences before and after alignment. In such cases the alignments will be generated within less than a second. The implementations are mathematically optimal but the libraries that we rely on do not scale well to longer sequences.

Use a specially designed software, that relies on heuristics, to perform large scale studies. These tools will operate orders of magnitude faster. Depending on your needs `blast`, `blat`, `mummer`, `minimap2`, `lastz`, `exonerate` will be far better suited for genome wide analyses.
 
## DNA alignment

Align the DNA corresponding to protein `S`

```{bash, comment=NA}
bio ncov:S ratg13:S --end 60 --align
```

## DNA alignment with 1 letter amino acid codes

```{bash, comment=NA}
bio ratg13:S ncov:S  --end 60 --align -1
```

Reading frame will follow the slice!

## DNA alignment with 3 letter amino acid codes

```{bash, comment=NA}
bio ratg13:S ncov:S  --end 60 --align -3
```

Reading frame will follow the slice!

## DNA alignment, tabular output

```{bash, comment=NA}
bio ncov:S ratg13:S --end 90 --align --table
```

## Align the translated regions

```{bash, comment=NA}
bio ncov:S ratg13:S --end 90 --translate --align 
```

## Align the protein corresponding to gene S

The protein sequence is fetched from the data (if exists) and is not a translated DNA. 

```{bash, comment=NA}
bio ncov:S ratg13:S --end 30 --protein --align 
```

The slice now applies to the protein sequence.

## Default alignment is global

With the default global alignment end gaps are have no penalty.

```{bash, comment=NA}
bio THISLINE ISALIGNED  -i --align
```

There is a strict mode that applies end gap penalties.

## Tabular output

All alignment may be formatted with tabular output

```{bash, comment=NA}
bio THISLINE ISALIGNED  -i --align --table
```

## Local alignment

Will produce all local alignments.

```{bash, comment=NA}
bio THISLINE ISALIGNED -i --align --local
```

## Global alignment

```{bash, comment=NA}
bio THISLINE ISALIGNED -i --align --global
```

## Semiglobal alignment

Same as zero endgap global but reports only the aligned region:

```{bash, comment=NA}
bio  THISLINE ISALIGNED -i --align --semiglobal
```

## Strict global alignment

Applies  end gap penalities.

```{bash, comment=NA}
bio THISLINE ISALIGNED -i --align --global --strict
```

