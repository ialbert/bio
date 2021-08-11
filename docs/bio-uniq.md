# `uniq.py`: find unique elements

The need to find unique elements within columns of different files is very common.

Thus when you install the `bio` package another script called `uniq.py` is also installed.
This software prints the unique elements from a column.

## Using `comm.py`

 If file 1 contains:

    A
    B
    C
    A
    B

    
then the command:

    uniq.py file_1.txt

will print:

    A
    B
    C

The flag `-c` used as:

    uniq.py -c file_1.txt
    
will print:

    2           A
    2           B
    1           C


`uniq.py` can be used from standard input:

    cat file_1.txt |  uniq.py -c

## Why does `uniq.py` exist?

We could use the UNIX construct:

    sort | uniq -c | sort -rn

the problem with the above is that the columns it prints are not tab separated. We may also use the entrez direct tool called:

    sort-uniq-count-rank

but for that `entrez-direct` must be installed.

## Additional utility

`uniq.py` can read different columns of a file and the delimiter may be changed as well. Read the second columns of three comma separated files:

    uniq.py -c -d , -f 2  file1 file2 file3

I don't usually advocate rewriting UNIX tools, in this case, writing a better `uniq` makes a lot of sense.

## Usage

```{bash, comment=NA}
uniq.py -h
```
