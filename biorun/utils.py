"""
Utilites funcions.
"""
import glob
import gzip
import json
import logging
import os
import re
import requests
import shutil
import sys
import tempfile
import time
from functools import wraps
from itertools import count, islice
from os.path import expanduser
from pprint import pprint

from biorun.libs.sqlitedict import SqliteDict

# The path to the current file.
__CURR_DIR = os.path.dirname(__file__)

# The default path to templates.
__TMPL_DIR = os.path.join(__CURR_DIR, "templates")

DATADIR = os.path.join(expanduser("~"), ".bio")

# Create the cache directory
os.makedirs(DATADIR, exist_ok=True)


def read_lines(stream, index=0):
    lines = filter(lambda x: not x.startswith("#"), stream)
    lines = map(lambda x: x.strip().split()[index], lines)
    lines = filter(None, lines)
    lines = list(lines)
    return lines


def parse_alias(fname):
    """
    Aliases hard to read accessions to easier to read names:
    NC_045512  wuhan-hu-1
    MN996532   raTG13
    """
    if not fname or not os.path.isfile(fname):
        return dict()

    stream = open(fname)
    lines = map(lambda x: x.strip(), stream)
    lines = filter(lambda x: not x.startswith("#"), lines)
    lines = filter(None, lines)
    lines = map(lambda x: x.split(), lines)
    lines = filter(lambda x: len(x) > 1, lines)
    pairs = map(lambda x: (x[0].strip(), x[1].strip()), lines)
    remap = dict(pairs)
    return remap


def is_int(text):
    try:
        int(text)
        return True
    except ValueError as exc:
        return False


def trim(text, size=3):
    """
    Trims a sequence to a length that is the largest multiple of size.
    """
    div, mod = divmod(len(text), size)
    subs = text[:div * size]
    return subs


def lower_case_keys(adict):
    return dict((k.lower(), v) for (k, v) in adict.items())


def safe_int_zero(text):
    try:
        return int(text)
    except ValueError as exc:
        return 0


def progress_bar(frac, barlen=30, null=' ', marker='=', head=">"):
    pos = int(frac * barlen)
    bar = marker * pos + head + null * int(barlen - pos)
    return bar


def open_db(table, fname, flag='c'):
    """
    Opens a connection to a data table.
    """

    conn = SqliteDict(fname, tablename=table, flag=flag, encode=json.dumps, decode=json.loads)
    return conn


CHUNK = 25000


def save_table(name, obj, fname, flg='w'):
    size = len(obj)
    table = open_db(table=name, fname=fname, flag=flg)
    for index, (key, value) in enumerate(obj.items()):
        table[key] = value
        if index % CHUNK == 0:
            perc = round(index / size * 100)
            print(f"*** saving {name} with {size:,} elements ({perc:.0f}%)", end="\r")
            table.commit()
    print(f"*** saved {name} with {size:,} elements (100%)", end="\r")
    print("")
    table.commit()
    table.close()


def plural(target, val=0, end='ies'):
    """
    Make the target string plural
    """

    output = target if val == 1 else f"{target[:-1]}{end}"

    return output


def response(url, params={}):
    # Open request to file
    r = requests.get(url, stream=True, params=params)
    try:
        # Check valid response status.
        r.raise_for_status()
    except Exception as exc:
        error(f"{exc}")
    return r


def download(url, dest_name, cache=False, params={}):
    """
    Downloads a URL into a destination
    """
    logger.info(f"downloading: {url}")

    # The file destination.
    if cache:
        path = os.path.join(DATADIR, dest_name)
    else:
        path = dest_name

    # Keep track of total size.
    total = 0

    # Open request
    url = url.replace('ftp:', 'http:') if url.startswith('ftp:') else url

    r = response(url=url, params=params)
    # Attempt to determine the download size.
    headers = lower_case_keys(r.headers)
    size = headers.get("content-length", 0)

    size = safe_int_zero(size)

    # How much data to process at a time.
    chunk_size = 1 * 1024 * 1024

    # The name of the file that will be stored
    fname = os.path.split(dest_name)[-1]

    # Create file only if download completes successfully.
    with tempfile.NamedTemporaryFile() as fp:

        # Iterate over the content and write to temp file
        for chunk in r.iter_content(chunk_size=chunk_size):
            total += len(chunk)

            if size:
                frac = 1 if total >= size else total / size
                perc = frac * 100
                bar = progress_bar(frac)
                print(f"*** downloading [{bar}] {fname} {human_size(size)} ({perc:.1f}%)", end="\r")
            else:
                print(f"*** downloading {fname} ({human_size(total)})     ", end="\r")

            fp.write(chunk)

        print("")

        # File creation completed.
        fp.seek(0)

        # Copy file to destination.
        shutil.copyfile(fp.name, path)

        # Progress notification.
        logger.info(f"saved to: {dest_name}")



def no_dash(alist):
    """

    """
    elems = list(filter(lambda x: x.startswith('-'), alist))
    if elems:
        msg = f"Invalid accessions: {elems}"
        error(msg)


def zero_based(start, end):
    """
    Shift to zero based coordinate system.
    """

    # Don't allow zero for start.
    if start == 0:
        error(f"start={start} may not be  zero")

    # Default value for start
    start = int(start) if start else 1

    # Default value for end.
    end = None if not end else int(end)

    # Shift the coordinates
    start = start - 1 if start > 0 else start
    return start, end


def human_size(num):
    for unit in ['B', 'KB', 'MB', 'GB']:
        if abs(num) < 1024.0:
            return "%.0f %s" % (num, unit)
        num /= 1024.0
    return "%.1f%s" % (num, '??')


def safe_int(text):
    try:
        return int(text)
    except ValueError as exc:
        logger.error(f"not a valid integer value: {text}")
        sys.exit()


def parse_number(text):
    """
    Parses a number from alternative representations: 100000, 100,000 or 100Kb or 100k all have the same representation.
    """
    if text in ('', None):
        return None

    text = str(text)
    text = text.lower()

    # Get rid of commas
    text = text.replace(",", '')

    if text.endswith("k") or text.endswith("kb"):
        value = safe_int(text.split("k")[0])
        text = f"{value * 1000}"

    if text.endswith("m") or text.endswith("mb"):
        value = safe_int(text.split("m")[0])
        text = f"{value * 1000 * 1000}"

    return safe_int(text)


def gz_write(fname, flag='wt'):
    """
    Shortcut to opening gzipped or regular files
    """
    stream = gzip.open(fname, flag, compresslevel=3) if fname.endswith(".gz") else open(fname, flag)
    return stream


def gz_read(fname, flag='rt'):
    """
    Shortcut to opening gzipped or regular files
    """
    stream = gzip.open(fname, flag) if fname.endswith(".gz") else open(fname, flag)
    return stream


def get_logger(name="bio", hnd=None, fmt=None, terminator='\n'):
    """
    Initializes a logger with a handler and formatter.
    """
    # Get the logger name.
    log = logging.getLogger(name)

    # Default logging level.
    log.setLevel(logging.WARNING)

    # The log handler.
    hnd = hnd or logging.StreamHandler()

    hnd.terminator = terminator

    # The logging formatter.
    fmt = fmt or logging.Formatter('*** %(message)s')

    # Add formatter to handler
    hnd.setFormatter(fmt)

    # Add handler to logger
    log.addHandler(hnd)

    return log


def set_verbosity(logger, level=1):
    """
    Sets the verbosity of the logger.
    """
    level = logging.DEBUG if level > 0 else logging.WARNING
    logger.setLevel(level)


# Initialize the logger.
logger = get_logger("main")


def error(msg, logger=logger, stop=True):
    """
    The default error handler
    """
    logger.error(f"{msg}")
    if stop:
        sys.exit(1)
