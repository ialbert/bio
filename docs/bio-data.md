# Data management {#bio-data}

The `bio` package helps maintain clarity and order among diverse datasets. When you obtain data with `bio` it becomes universally available to tool in `bio` regardless of your location on the filesystem.

## Getting data from NCBI (--fetch)

The `--fetch` command downloads data identified via accession numbers from NCBI then stores 
this data in a storage directory (`~/.bio`). All subsequent commands in the `bio` package can seamlessly access the stored  data from any location and would not need to connect to the internet to use it.

    # Run fetch in verbose mode.
    bio NC_045512 --fetch -v
    
Running the fetch command the subsequent times for the same accession number will not connect to the internet again, it will use the existing data instead. Use the `--fetch --update` (see later) to force a re-downloading of data from NCBI. 

Most commands can operate on multiple accession numbers at a time.

    bio NC_045512 MN996532 --fetch
    
Depeding on the datasize and internet connection speed the first fetch may take from seconds (1 second SARS-COV-2, 675 lines) to minutes (7 minutes for Human Chromosome 1, 4.9 million lines) or possibly longer.

The internal, gzip compressed, JSON based representation used by `bio` is simple and efficient. The 330MB GenBank file of chromosome 1 (human genome) takes just 67MB to store in our representation. More importantly `bio` can read and convert the stored information to a `fasta` format containing 253 million basepairs in just 6 seconds:

    time bio NC_000001 --fasta | wc -c
    253105766
    
    real    0m6.189s

For shorter genomes the conversion times are proportionally shorter.

## Rename  (--rename)

Accession numbers are tedious to handle. Almost always we rename data to be meaningful.

    bio NC_045512 --fetch --rename ncov

the command above will store the data under the name `ncov`. Within the data the sequence will be labeled with its version number `NC_045512.2`. You may change both the name and sequence id:

    bio NC_045512 --fetch --rename ncov --seqid ncov

##  Listing the storage (--list)

Each time data is fetched is will be stored locally. The `bio` package will look in this storage no matter what directory you run it from. To list the content of the storage write:

    bio --list

it prints:

    22K   ncov     Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    19K   ratg13   Bat coronavirus RaTG13, complete genome

## Delete data (--delete)

To drop data from storage use:

    bio ncov --delete
    
This command will only drop the JSON representation not the downloaded GenBank if exists.
If you want to update the original data use the `--update` parameter.

## Update data (--update)   
    
To force fetch to download data that already seems to be present  do:

    bio NC_045512 --fetch --update

Note that you can't update a renamed sequence. At that point the original accession number is lost. You can however fetch, update and rename all in one go like so:

    bio NC_045512 --fetch --update --rename ncov --seqid ncov

There is a predetermined order of operations, thus it does not matter in what order you parameters. For example `--delete` would take place first before the `--fetch` and so on.

## View data

The default action is to view the stored data.  Locally the data is stored in a JSON format that makes processing it much faster than the original GenBank yet has no loss of information:

    bio ncov | head 
         
## Reading files

`bio` can also read and process data from local files. The software will attempt to detect the type of the file from the file extention. .gb, .gbk or .genebank for GenBank, .fa, .fasta for FASTA. Here is how to turn a genbank file into gff:

    bio  mydata.gb --gff 
 
## Data from command line

In addition, there is a so called *interactive* input of data (`-i`) where the data can be listed at the command line:

    bio ATGAATATATAC -i --translate
   
The above command will operate on the sequence as if it were stored in a FASTA file, the above prints:
    
    >1 translated
    MNIY

Interative more works for exploration/demonstration. Suppose you wanted to see how the same DNA sequence would be translated to different peptides in different reading frames:

```{bash, comment=NA}
bio ATGAATATATACT -i --translate --start 1
```

versus

```{bash, comment=NA}
bio ATGAATATATACT -i --translate --start 2
```

or you can explore alignments:

```{bash, comment=NA}
    bio THISLINE ISALIGNED --align -i 
```


