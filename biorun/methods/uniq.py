import codecs
import csv
import sys
from collections import defaultdict

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

@plac.opt('field', "field index (1 by default)", type=int)
@plac.opt('delim', "delimiter (guess by default)")
def main(field=1, delim=''):
    # Input stream.
    stream = sys.stdin

    delim = decode(delim)

    if delim:
        reader = csv.reader(stream, delimiter=delim)
    else:
        # Figure out the dialect from the file.
        try:
            dialect = csv.Sniffer().sniff(stream.read(1024))

            # Rewind stream to the start.
            stream.seek(0)

            # Read the stream.
            reader = csv.reader(stream, dialect)

        except Exception as exc:
            reader = csv.reader(stream, delimiter="\t")

    # Stores the counter.
    store = defaultdict(int)

    # Fill the
    idx = field - 1
    for row in reader:

        try:
            key = row[idx]
        except IndexError as exc:
            key = ''

        store[key] += 1

    pairs = [(v, k) for k, v in store.items()]

    pairs.sort(key=lambda x: (-x[0], x[1]))
    for v, k in pairs:
        print(f"{v}\t{k}")

@nointerrupt
def run():
    """
    Entry point for the script.
    """
    plac.call(main)


if __name__ == '__main__':
    run()
