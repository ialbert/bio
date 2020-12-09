# comm.py

When you install the `bio` package you get another script called `comm.py` is also installed.
It is a tool that prints the common elements from two files. If file 1 contains:

    A
    B
    C

and file 2 contains 

    A
    C
    D
    
the `comm.py file1 file2` will print:

    A
    C
    
That's it. These are the elements in common in the first column of both files.    

## Other features

`comm.py` has a number of convenience features, it can:

1. print elements unique to file 1: `B`
1. print elements unique to file 2: `D`
1. print the union of elements" `A`, `C`, `B`, `D`
1. read different columns of the files 
1. read CSV and tab-delimited files

The content for either file may come from standard input. In that case the `-` symbol should be used instead of file name.

## Rationale

The need to find identical elements within columns of different files is surprisingly common (pun intended).

We could use the UNIX tool called `comm` to find common or distinct elements. When used properly `comm` allows you to answer a wide variety of interesting questions. Unfortunately using `comm` properly is no easy task. First the values must be on a single column and must be sorted. Then the `comm` command may feel exceedingly counter-intuitive. Instead of telling it what we want, we have to tell it what we don’t want (what columns to suppress). It is completely backwards of how I like to think. While I don't usually advocate rewriting UNIX tools, in this case, writing a better `comm` makes a lot of sense.

## Limitations

With `comm.py` most operations will be quicker to do, simpler to perform and easier to understand. The primary limitation of `comm.py` vs `comm` is that `comm.py` loads all elements into memory. Once the number of elements passes about 1 million `comm.py` will be noticeably and increasingly slower than `comm`. Under 1 million items using `comm.py` will work fine. For most usecases `comm.py` works exceedingly well.

## Usage

```{bash, comment=NA}
comm.py -h
```
