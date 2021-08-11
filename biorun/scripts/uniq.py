import codecs
import csv
import sys
from collections import defaultdict
from biorun import utils
import biorun.libs.placlib as plac

def nointerrupt(func):
    """
    Intercept keyboard interrupts.
    """
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except KeyboardInterrupt:
            sys.exit(0)

    return wrapper

def decode(text):
    """
    Recognize string encodings: \t etc
    """
    return codecs.decode(text, 'unicode_escape')

@plac.flg('count', "produce counts")
@plac.opt('field', "field index (1 by default)", type=int)
@plac.opt('delim', "delimiter (tab by default)")
@plac.pos('fnames', "file names")
def main(field=1, delim='', count=False, *fnames):

    streams = utils.get_streams(fnames)

    # Find the delimiter.
    delim = decode(delim) if delim else "\t"

    # Initialize the counter.
    store = defaultdict(int)

    for stream in streams:

        # Parse the file.
        reader = csv.reader(stream, delimiter=delim)

        # Fill the dictionary.
        idx = field - 1
        for row in reader:
            try:
                key = row[idx]
            except IndexError as exc:
                key = ''

            store[key] += 1


    if count:
        # Produce sorted counts.
        pairs = [(v, k) for k, v in store.items()]
        pairs.sort(key=lambda x: (-x[0], x[1]))
        for v, k in pairs:
            print(f"{v}\t{k}")
    else:
        for k in store:
            print(f"{k}")

@nointerrupt
def run():
    """
    Entry point for the script.
    """
    plac.call(main)


if __name__ == '__main__':
    run()
