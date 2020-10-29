# bio: data

The `bio` package solves the ongoing struggle of how to maintain sanity among diverse datasets.

When you obtain data with `bio` it becomes universally available to all tools in the package.

### Getting data (--fetch)

The `--fetch` command downloads data identified via accession numbers from NCBI then stores 
this data in a cache directory. All subsequent commands in the `bio` package can seamlessly access the cached 
data from any location and would not need to connect to the internet.

    # Get data for a single accession number.
    bio NC_045512 --fetch
    
    # Get data for multiple accession numbers.
    bio NC_045512 MN996532 --fetch
    
### Rename  (--rename)

Accession numbers are tedious to handle. Almost always we rename data to be meaningful

    bio NC_045512 --fetch --rename SARS2

the command above will rename store the data under the name `SARS2`. Within the data the sequence will still be labeled as `NC_045512`. We can change the internal id of the sequence with:

    bio NC_045512 --fetch --rename SARS2 --seqid SARS2

Now, not just the file is called `SARS2` but the sequence id is also set to `SARS2` as well.

###  Listing the storage (--list)

Each time data is fetched is stored in the cache system. The `bio` package will always find it no matter what directory you run it from. To list the content of the storage write:

    bio --list

it prints:

    22K  SARS2   Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome

### Delete data (--delete)

To drop data from storage use:

    bio SARS2 --delete
    
### Update data (--delete --fetch)   
    
When data is updated a more recent version may be available. Delete and refetch to update data:

    bio NC_045512 --delete --fetch

There is a builtin order of operations, even if you wrote it backwards (`--fetch --delete`) the deletion would take place first. You can chain all other actions as well:

    bio NC_045512 --delete --fetch --rename SARS2 --seqid SARS2

### View data

The default action is to view the stored data.  When the data is stored it is transformed into a format that makes processing much faster than the original GenBank yet has no loss of information. The new format is so called JSON:

    bio NC_045512 | head 
 
View data that is tagged with the `S` gene:

    bio NC_045512 | head
    
## Example output

```{bash, comment=NA}
bio NC_045512 | head
```

