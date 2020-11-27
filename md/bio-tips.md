# General tips {#tips}

Effort has been devoted to make the command line a bit more user friendly.

## Interactive mode

Passing the `-i` flag allows data to be passed from command line. For example:

```{bash, comment=NA}
bio --translate -i ATGATTATATATA 
```

Note how the input was read as parameters from the  command line. We make use of this feature when exploring simple data in an explicit way.

## The paramter format is forgiving

You may use single or double dashes on parameters:

    bio ncov --fasta --end 100
    
or:

    bio ncov -fasta -end 100
    
both will work the same way.

## The coordinate system is 1 based

Coordinates are 1 based (inclusive on both ends) identical to GFF coordinate formats.

Numbers for start and end coordinates may be written in human friendly form, like so: `5000` or `5,000` or `5K` or `5KB`.
  
## The order of operations is pre-determined

You may combine multiple parameters, in that case each condition will be applied sequentially in a internally detrermined order
that is independent of the order the parameters are listed. 

For example when using a `--start` and `--end` and `--translate` the selection by start and end takes place on the DNA then the resulting sequence is translated into aminoacids.
 
The same start, end combo followed by `--protein` applies the slice on the protein sequences as aminoacids.
  
## You may use multiple accession numbers
   
Many commands allow using multiple accession numbers, in that case the operations will take place sequentially on each.

    bio NC_045512 MN996532 --fetch 
  
## Generate more verbose outputs

Use the `-v` flag to produce verbose outputs for each command.
