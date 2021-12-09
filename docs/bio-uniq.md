# bio uniq: find unique elements {#bio-uniq}

The need to find unique elements within columns of different files is very common.


## Using `comm.py`

 If file 1 contains:

    A
    B
    C
    A
    B

    
then the command:

    bio uniq file_1.txt

will print:

    A
    B
    C

The flag `-c` used as:

    bio uniq  -c file_1.txt
    
will print:

    2           A
    2           B
    1           C


`uniq.py` can be used from standard input:

    cat file_1.txt |  bio uniq -c

## Why does `uniq.py` exist?

We could use the UNIX construct:

    sort | uniq -c | sort -rn

the problem with the above is that the columns it prints are not tab separated. We may also use the entrez direct tool called:

    sort-uniq-count-rank

but for that `entrez-direct` must be installed.


In addition `bio uniq` can read different columns of a file plus the delimiter may be changed as well. To find the unique elements listed in the second column of three comma separated files:

    bio uniq -c -d , -f 2  file1 file2 file3

I don't usually advocate rewriting UNIX tools, in this case, writing a better `uniq` makes a lot of sense.

## Usage

```{bash, comment=NA}
bio uniq -h
```
