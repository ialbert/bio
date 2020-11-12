"""
Utilites funcions.
"""
import sys, os, re, tempfile, gzip, glob, shutil

from itertools import count, islice
from functools import wraps
import time
import logging
from os.path import expanduser

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


def save_stream(stream, fname, trigger=10000, file=sys.stdout, flag='wt'):
    """
    Write a input 'stream' as the fname filename
    """
    # Save a stream into file and print progess.
    # Use a temporary file in case the process fails.
    # We don't want a failed cache file.

    tmp = tempfile.NamedTemporaryFile(mode="w+t")

    sequence = count(1)
    for index, line in zip(sequence, stream):
        if (index % trigger) == 0:
            print(f"*** downloaded {index:,d} lines\r", file=file, end='')
        tmp.write(line)

    print()
    # Not sure if this is needed. Can't hurt.
    tmp.flush()

    # Rewind temporary file to beginning
    tmp.seek(0)

    # Copy over the content from the temporary file to the final gzipped file destination.
    out = gzip.open(fname, flag) if fname.endswith(".gz") else open(fname, flag)
    for line in tmp:
        out.write(line)
    out.close()
    tmp.close()

    logger.info(f"saved {fname}")

    return


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


def symlink(src, dest):
    """
    Creates a symlink.
    """

    # Allow link replacement.
    if os.path.islink(dest):
        os.remove(dest)

    # Don't link to files.
    if os.path.isfile(dest):
        logger.error(f"invalid link destination {dest}")

    os.symlink(src, dest)

# Initialize the logger.
logger = get_logger("main")

def error(msg, logger=logger):
    """
    The default error handler
    """
    logger.error(f"error: {msg}")
    sys.exit(1)
