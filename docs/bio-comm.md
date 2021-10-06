# bio comm: find common values {#bio-comm}

The need to find identical elements within columns of different files is surprisingly common (pun intended).

Thus `bio` packages another utility called `comm`. It is a tool that prints the common elements from two files.'

## Using `bio comm`

 If file 1 contains:

    A
    B
    C

and file 2 contains 

    A
    C
    D
    
then the command:

    bio comm file1 file2

will print the common elements present in the first column of both files:

    A
    C

This is the main usecase of the `bio comm` software.

## Other features

`bio comm` has a number of convenience parameters:

* `-1` will print elements unique to file 1: `B`
* `-2` will print elements unique to file 2: `D`
* `-3` will print the union of elements: `A`, `C`, `B`, `D`
* `-x 1` reads a different column from file 1
* `-y 1` reads a different column from file 2
* `-t` treats the files as tab delimited rather than CSV

The content for either file may come from standard input. In that case the `-` symbol should be used instead of file name.

## Why does `bio comm` exist?

We could use the UNIX tool called `comm` to find common or distinct elements. When used properly `comm` allows you to answer a wide variety of interesting questions.

Unfortunately using `comm` properly is no easy task.

First for `comm` to work the values must be on a single column and must be sorted. Then instead of telling `comm` what we want, we have to tell it what we don’t want (what columns to suppress). That usage is completely backwards of how I like to think.

I don't usually advocate rewriting UNIX tools, in this case, writing a better `comm` makes a lot of sense.

## Potential limitations

With `bio comm` most operations will be quicker to do, simpler to perform and easier to understand. The primary limitation of `bio comm` vs `comm` is that `bio comm` loads all elements into memory.

For most use-cases `bio comm` will work exceedingly well.

## Usage

```{bash, comment=NA}
bio comm -h
```
