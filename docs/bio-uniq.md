# uniq.py

The need to find unique elements within columns of different files also very common.

Thus when you install the `bio` package another script called `uniq.py` is also installed.
It is a tool that prints the unique elements from a column.

## Using `comm.py`

 If file 1 contains:

    A
    B
    C
    A
    B

    
then the command:

    cat foo | uniq.py

will print:

    A
    B
    C

The flag `-c` used as:

    cat foo | uniq.py -c
    
will print:

    2           A
    2           B
    1           C

## Why does `uniq.py` exist?

We could use the UNIX construct

    sort | uniq -c | sort -rn

or the entrez direct tool called:

    sort-uniq-count-rank

I don't usually advocate rewriting UNIX tools, in this case, writing a better `uniq` makes a lot of sense.

## Usage

```{bash, comment=NA}
uniq.py -h
```
