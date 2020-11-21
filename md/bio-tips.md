# General tips {#tips}

## Interactive mode

Passing the `-i` flag allows data to be passed from command line. For example:

```{bash, comment=NA}
bio --translate -i ATGATTATATATA 
```

Note how the input was read as parameters from the  command line. We make use of this feature when exploring simple data in an explicit way.
Here we're using a different phase to translate the same sequence above.

```{bash, comment=NA}
bio --translate -i ATGATTATATATA --start 2
```

## The command line is more forgiving:

You may write short or long form commands:

    bio ncov --fasta --end 100
    
or you may write:

    bio ncov -fasta -end 100
    
both will work the same way.

## The coordinate system is 1 based

Coordinates are 1 based (inclusive on both ends) identical to GFF coordinate formats.

## The order of operations is pre-determined

You may combine multiple parameters, in that case each condition will be applied sequentially in a internally detrermined order
that is independent of the order the parameters are listed. 

For example when using a `--start` and `--end` and `--translate` the selection by start and end takes place on the DNA then the resulting sequence is translated into aminoacids.
 
The same start, end combo followed by `--protein` applies the slice on the protein sequences as aminoacids.
  
## You may use multiple accession numbers
   
Many commands allow using multiple accession numbers, in that cose the operations will take place sequentially on each.

## Generate more verbose outputs

Use the `-v` flag to produce verbose outputs for each command.
