"""
Utilites funcions.
"""
import sys, os, re, tempfile, gzip, glob, shutil

from itertools import count, islice
from functools import wraps
import time
import logging
from os.path import expanduser
from biorun import const

# The path to the current file.
__CURR_DIR = os.path.dirname(__file__)

# The default path to templates.
__TMPL_DIR = os.path.join(__CURR_DIR, "templates")

DATADIR = os.path.join(expanduser("~"), ".bio")

# Create the cache directory
os.makedirs(DATADIR, exist_ok=True)


def time_it(func):
    @wraps(func)
    def timer(*args, **kwargs):
        units = "seconds"
        start = time.time()
        try:
            return func(*args, **kwargs)
        finally:
            end = time.time()
            diff = int(round((end - start), 1)) or 0.1
            if diff > 120:
                diff, units = diff / 60, "minutes"
            logger.info(f"{func.__name__} runtime: {diff} {units}")

    return timer


def maybe_prot(name):
    """
    Name may be a proteing accession
    """
    return name[:2] in const.NCBI_PROTEIN_CODES


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
    subs = text[:div*size]
    return subs

def maybe_ncbi(text):
    """
    Guesses that a text is a valid NCBI accession
    """
    code = text[:2]

    # Upper case letters
    upper = code == code.upper()

    # Get rid of underscores and version number
    rest = text[2:].replace("_", "").split(".")[0]

    # Sizes between 5 and 9
    size = 5 <= len(rest) <= 9

    return upper and size and is_int(rest)


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


def sizeof_fmt(num, suffix=''):
    for unit in ['', 'K', 'M', 'G']:
        if abs(num) < 1024.0:
            return "%.0f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, '??', suffix)


def safe_int(text):
    try:
        return int(text)
    except ValueError as exc:
        error(f"not an integer value: {text}")


def parse_number(text):
    """
    Parses a number from alternative representations: 100000, 100,000 or 100Kb or 100k all have the same representation.
    """
    text = text.lower()

    text = text.replace(",", '')

    if text.endswith("k") or text.endswith("kb"):
        value = safe_int(text.split("k")[0])
        text = f"{value * 1000}"

    if text.endswith("m") or text.endswith("mb"):
        value = safe_int(text.split("m")[0])
        text = f"{value * 1000 * 1000}"

    return text


def save_stream(stream, fname, trigger=10000, file=sys.stderr, flag='wt'):
    """
    Write a input 'stream' as the fname filename
    """
    # Save a stream into file and print progess.
    # Use a temporary file in case the process fails.
    # We don't want a failed cache file.

    tmp = tempfile.NamedTemporaryFile(mode="w+t")

    index = 0

    sequence = count(1)

    for index, line in zip(sequence, stream):
        if (index % trigger) == 0:
            print(f"*** downloaded {index:,d} lines\r", file=file, end='')
        tmp.write(line)

    print(f"*** downloaded  {index:,d} lines", file=file)

    # Not sure if this is needed. Can't hurt.
    tmp.flush()

    # Rewind temporary file to beginning
    tmp.seek(0)

    # Copy over the content from the temporary file to the final gzipped file destination.
    out = gz_write(fname)
    for line in tmp:
        out.write(line)
    out.close()
    tmp.close()

    logger.info(f"saved {fname}")

    return


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


def get_logger(name, hnd=None, fmt=None, terminator='\n'):
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


def symlink(src, dst):
    """
    Creates a symlink.
    """

    # Allow link replacement.
    if os.path.islink(dst):
        os.remove(dst)

    # Don't link to files.
    if os.path.isfile(dst):
        logger.error(f"invalid link destination {dst}")

    os.symlink(src, dst)


# Initialize the logger.
logger = get_logger("main")


def error(msg, logger=logger):
    """
    The default error handler
    """
    logger.error(f"error: {msg}")
    sys.exit(1)
