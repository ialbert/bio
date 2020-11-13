# Sequence alignments {#align}


# Align the extracted protein.
bio align ncov:S ratg13:S --end 80 > align-dna-s.txt

# Align the extracted protein.
bio align ncov:S ratg13:S --protein > align-protein-s.txt

# Align the translated regions.
bio align ncov:S ratg13:S --end 80 --translate > align-translated-s.txt

# Test alignments
bio align THISLINE ISALIGNED  -i > align-local.txt

# Default alignment is local.
bio align THISLINE ISALIGNED  -i > align-local.txt

# Global alignment.
bio align THISLINE ISALIGNED -i --mode global > align-global.txt

# Semiglobal alignment.
bio align THISLINE ISALIGNED -i --mode semiglobal > align-semiglobal.txt

Strict global.

```{bash, comment=NA}
bio align THISLINE ISALIGNED -i --mode strictglobal
```


```{bash, comment=NA}
bio align -i THISLINE ISALIGNED
```
