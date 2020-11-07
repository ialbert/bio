# bio: data {#data}

The `bio` package solves the ongoing struggle of how to maintain sanity among diverse datasets.

When you obtain data with `bio` it becomes universally available to all tools in the package.

### Local data

There is an  automated data storage in `bio` that makes using data from NCBI quite convenient. But before we even go ther let's make it clear that `bio` can read and process data from local file just as well. If you have a genbank or fasta file you can use that as input:

    bio  mydata.gb --fasta 

Moreover there is a so called command line input of data (`-i`) where the data is read as listed on the command line:

    bio ATGAATATATAC -i --translate
   
will operate on the sequence listed as if it were a DNA sequence:     
    
    >S1 translated DNA
    MNIY

### Getting data from NCBI (--fetch)

The `--fetch` command downloads data identified via accession numbers from NCBI then stores 
this data in a storage directory (`~/.bio`). All subsequent commands in the `bio` package can seamlessly access the stored  data from any location and would not need to connect to the internet to use it.

    # Get data for a single accession number.
    bio NC_045512 --fetch
    
Running the fetch command once you already have the data will not connect to the internet again.

Most commands you can operate on multiple accession numbers at a time.

    bio NC_045512 MN996532 --fetch
    
There will be commands like `--rename` where it makes no sense to apply the operation on multiple data at the same time. In those cases only the first accession number is acted upon.

### Rename  (--rename)

Accession numbers are tedious to handle. Almost always we rename data to be meaningful.

    bio NC_045512 --fetch --rename ncov

the command above will store the data under the name `ncov`. Within the data the sequence will still be labeled as `NC_045512`. You may change both the name and sequence id:

    bio NC_045512 --fetch --rename ncov --seqid ncov
    bio NC_045512 --fetch --rename ratg13 --seqid ratg13
    
Now, not only is your data called `ncov` but the sequence id inside the data is also set to `ncov`.

###  Listing the storage (--list)

Each time data is fetched is will be stored locally. The `bio` package will look in this storage no matter what directory you run it from. To list the content of the storage write:

    bio --list

it prints:

    22K   ncov     Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    19K   ratg13   Bat coronavirus RaTG13, complete genome

### Delete data (--delete)

To drop data from storage use:

    bio ncov --delete
    
This command will only drop the JSON representation not the downloaded GenBank if exists.
If you want to update the original data use the `--update` parameter.

### Update data (--update)   
    
To force fetch to download again the data that you already have:

    bio NC_045512 --fetch --update

Note that you can't update a renamed sequence. At that point the original accession number is lost. You can however fetch, update and rename all in one go like so:

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


