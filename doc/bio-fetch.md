# bio: data {#data}

The `bio` package solves the ongoing struggle of how to maintain sanity among diverse datasets.

When you obtain data with `bio` it becomes universally available to all tools in the package.

### Getting data (--fetch)

The `--fetch` command downloads data identified via accession numbers from NCBI then stores 
this data in a cache directory. All subsequent commands in the `bio` package can seamlessly access the cached 
data from any location and would not need to connect to the internet.

    # Get data for a single accession number.
    bio NC_045512 --fetch
    
For most commands you can operate on multiple accession numbers at a time.

    bio NC_045512 MN996532 --fetch
    
There will be commands like `--rename` where it makes no sense to apply the operation on multiple data at the same time. In those cases only the first accession number is acted upon.

### Rename  (--rename)

Accession numbers are tedious to handle. Almost always we rename data to be meaningful

    bio NC_045512 --fetch --rename ncov

the command above will rename store the data under the name `ncov`. Within the data the sequence will still be labeled as `NC_045512`. We can change the internal id of the sequence with:

    bio NC_045512 --fetch --rename ncov --seqid ncov
    bio NC_045512 --fetch --rename ratg13 --seqid ratg13
    
Now, not just the file is called `ncov` but the sequence id is also set to `ncov` as well.

###  Listing the storage (--list)

Each time data is fetched is stored in the cache system. The `bio` package will always find it no matter what directory you run it from. To list the content of the storage write:

    bio --list

it prints:

    22K   ncov     Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    19K   ratg13   Bat coronavirus RaTG13, complete genome

### Delete data (--delete)

To drop data from storage use:

    bio ncov --delete
    
This command will only drop the JSON representation not the downloaded GenBank if exists.
If you need to update the original data use the `--update` parameter.

### Update data (--update)   
    
To get the latest from NCBI update data:

    bio NC_045512 --fetch --update

Note that we can't update a renamed sequence, at that point the original accession number is lost. What we can do is fetch and update the accession number then rename all in one line like so:

    bio NC_045512 --fetch --update --rename ncov --seqid ncov

Since the update needs internet access, and depending on datasize waiting for download use it only if you have reason to believe the data has changed.

There is a builtin order of operations, does not matter what order you list commands. For example deletion would take place first before the fetch and so on.

### View data

The default action is to view the stored data.  Locally the data is stored in a JSON format that makes processing it much faster than the original GenBank yet has no loss of information:

    bio ncov | head 
 
You can subselect and view sections of the data, for example the `S` gene:

    bio ncov:S | head
    
Many other combinations are valid:

    bio ncov --type CDS 
    
## Example output

```{bash, comment=NA}
bio ncov | head
```


