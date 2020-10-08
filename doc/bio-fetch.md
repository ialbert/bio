# bio: fetch

The `fetch` command downloads data identified via accession numbers from NCBI then stores 
this data in a cache directory. All subsequent commands in the `bio` package can seamlessly access the cached 
data from any location and would not need to connect to the internet.

    # Get data for a single accession number.
    
   bio fetch NC_045512  
    # Get data for multiple accession numbers.
    bio fetch NC_045512 MN996532
 
You may then use `view` to investigate the cached data:

    bio view NC_045512 | head 
    
## File input
   
The accession numbers to `fetch` may also come from a tab delimited file where the 
content of the first column will interpreted as accession numbers.

```
#
# Commented lines will be ignored.
#
# bio fetch will read the first column of the file
#
AF086833
MN996532
```

The file above could be used as:

    bio fetch acc.txt | head

## Updating cached data

Pass the `-u` flag to have `fetch` obtain the latest version (if that exists) for the data under the accession number:

    bio fetch NC_045512 -u

## Example output

```{bash, comment=NA}
bio fetch NC_045512 | head
```

## Command line help
```{bash, comment=NA}
bio fetch -h
```
