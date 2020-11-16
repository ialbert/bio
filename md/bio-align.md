# Sequence alignments {#align}


```{bash, comment=NA}
# Align the extracted protein.
bio ncov:S ratg13:S --end 90 --align
```

# Align the translated regions.

```{bash, comment=NA}
bio ncov:S ratg13:S --end 90 --translate --align 
```

Default alignment is semiglobal.

```{bash, comment=NA}
bio THISLINE ISALIGNED  -i --align
```

Local alignment.

```{bash, comment=NA}
bio THISLINE ISALIGNED -i --align --local
```

Global alignment.

```{bash, comment=NA}
bio THISLINE ISALIGNED -i --align --global
```

Semiglobal alignment

```{bash, comment=NA}
bio  THISLINE ISALIGNED -i --align --semiglobal
```
