# bio: tips {#tips}

## More forgiving command line

You may write short or long form commands:

    bio ncov --fasta --end 100
    
or you may write:

    bio ncov -fasta -end 100
    
both will work the same way.

## One based coordinate system

Coordinates are 1 based (inclusive on both ends) identical to GFF coordinate formats.

## Order of operations

You may also combine multiple parameters, in that case each condition will be applied sequentially in a predermined order
that is independent of the order the parameters are listed. 

For example when using a `--start` and `--end` and `--translate` the selection by start and end takes place on the DNA then the resulting sequence is translated into aminoacids.
 
The same start, end combo followed by `--protein` applies the slice on the protein sequences as aminoacids.
  
## Multiple accession numbers
   
When using multiple accession numbers the operations will take place sequentially on each.

## More verbose outputs

Use the `-v` flag to produce verbose outputs for each command.
