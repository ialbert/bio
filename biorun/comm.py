"""
A replacement for the unix comm tool.

Differences:

- files do not need to be sorted
- data may be fetched from different columns
- produces common elements by default

"""
import csv, sys, os
import biorun.libs.placlib as plac
from pathlib import Path

FILL_VALUE = ''

UNIQ1, UNIQ2, ISECT, UNION = range(4)

def process(stream1, stream2, delimiter, idx1, idx2, show):
    """
    Processes the files and prints the output
    """

    def parse(stream, idx):
        """
        Returns the value of a column at the column index.
        """

        # Skip comment lines
        stream = filter(lambda x: not x.startswith('#'), stream)

        # Ignore empty lines.
        stream = filter(lambda x: x.strip(), stream)

        # Format the stream.
        stream = csv.reader(stream, delimiter=delimiter)

        # Generate empty values on missing columns.
        for row in stream:
            try:
                yield (row[idx], None)
            except IndexError as exc:
                yield ('', None)

    # Make dictionaries, will maintain original item order.
    store1 = dict(parse(stream1, idx=idx1))
    store2 = dict(parse(stream2, idx=idx2))

    # Generate the various groupings.
    isect = [key for key in store1.keys() if key in store2]
    uniq1 = [key for key in store1.keys() if key not in store2]
    uniq2 = [key for key in store2.keys() if key not in store1]
    union = isect + uniq1 + uniq2

    # Select output based on flags.
    if show == UNIQ1:
        stream = uniq1
    elif show == UNIQ2:
        stream = uniq2
    elif show == UNION:
        stream = union
    else:
        stream = isect

    # Print the output
    for line in stream:
        print(line)

def get_stream(fname):
    """
    Turns filename into a stream.
    """
    if fname == '-':
        return sys.stdin

    if not os.path.isfile(fname):
        print(f"file not found: {fname}")
        sys.exit(1)

    return open(fname)

@plac.pos('file1', "input file 1", type=Path)
@plac.pos("file2", "input file 2", type=Path)
@plac.flg('uniq1', "prints elements unique to file 1", abbrev="1")
@plac.flg('uniq2', "prints elements unique to file 2", abbrev="2")
@plac.flg('union', "prints elements present in both files", abbrev="3")
@plac.flg('tab', "tab delimited (default is csv)", abbrev="t")
@plac.opt('col1', "column index for file 1 [default=1]", abbrev="x", type=int)
@plac.opt('col2', "column index for file 2 [default=1]", abbrev="y", type=int)
def run(file1, file2, uniq1=False, uniq2=False, union=False, tab=False, col1=1, col2=1):
    """
    A better 'comm' command. Prints elements common from columns from two files.
    """
    delimiter = "\t" if tab else ","

    idx1 = col1 - 1
    idx2 = col2 - 1

    # Figure out what the mode of operation is.
    show = ISECT
    show = UNIQ1 if uniq1 else show
    show = UNIQ2 if uniq2 else show
    show = UNION if union else show

    if not os.path.isfile(file1):
        print(f"file not found: {file1}")
        sys.exit(1)

    # Get a stream for each file
    stream1 = get_stream(file1)
    stream2 = get_stream(file2)

    # Process the file.
    process(stream1=stream1, stream2=stream2, delimiter=delimiter, idx1=idx1, idx2=idx2, show=show)


def main():
    """
    Entry point for the script.
    """
    if len(sys.argv)<2:
        sys.argv.append("-h")
    plac.call(run)


if __name__ == '__main__':
    main()
