"""
Utilites funcions.
"""
import sys, os, tempfile, gzip, glob, shutil

from itertools import count, islice
from functools import wraps
import time
import jinja2
import logging
from os.path import expanduser

# The path to the current file.
__CURR_DIR = os.path.dirname(__file__)

# The default path to templates.
__TMPL_DIR = os.path.join(__CURR_DIR, "templates")

DATADIR = os.path.join(expanduser("~"), ".bio")

# Create the cache directory
os.makedirs(DATADIR, exist_ok=True)

GENBANK, FASTA, GFF, BED, SAM, BAM = "genbank", "fasta", "gff", "bed", "sam", "bam"

TYPE_BY_EXTENSION = {
    "gb": GENBANK,
    "gbk": GENBANK,
    "genbank": GENBANK,
    "fa": FASTA,
    "fasta": FASTA,
    "bed": BED,
    "gff": GFF,
    "sam": SAM,
    "bam": BAM,
}




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
                diff, units = diff/60, "minutes"
            logger.info(f"{func.__name__} runtime: {diff} {units}")

    return timer

def guess_type(path):
    """
    Attempts to guess a file type from an extension.
    """
    name, ext = os.path.splitext(path)
    ext = ext.lower()
    ftype = TYPE_BY_EXTENSION.get(ext, "")
    return ftype



def sizeof_fmt(num, suffix=''):
    for unit in ['', 'K', 'M', 'G']:
        if abs(num) < 1024.0:
            return "%.0f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, '??', suffix)

def save_stream(stream, fname, trigger=50000):
    """
    Write a input 'stream' as the fname filename.
    """
    # Save a stream into file and print progess.
    # Use a temporary file in case the process fails.
    # We don't want a failed cache file.
    tmp = tempfile.NamedTemporaryFile(mode="w+t")

    sequence = count(1)
    for index, line in zip(sequence, stream):
        if (index % trigger) == 0:
            logger.info(f"wrote {index:,d} lines")
        tmp.write(line)

    # Not sure if this is needed. Can't hurt.
    tmp.flush()

    # Rewind temporary file to beginning
    tmp.seek(0)

    # Copy over the content from the temporary file to the final gzipped file destination.
    out = gzip.open(fname, 'wt')  if fname.endswith(".gz") else open('wt')
    for line in tmp:
        out.write(line)
    out.close()
    tmp.close()

    # User friendly file size.
    #fsize = sizeof_fmt(os.path.getsize(fname))
    #logger.debug(f"saved {fsize} data to {fname}")
    return

def get_logger(name, hnd=None, fmt=None, terminator='\n'):
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
    level = logging.DEBUG if level > 0 else logging.WARNING
    logger.setLevel(level)

def symlink(src, dest):
    """
    Creates a symlink.
    """

    if os.path.islink(dest):
        os.remove(dest)

    if os.path.isfile(dest):
        logger.error(f"invalid link destination {dest}")

    os.symlink(src, dest)

# Initialize the logger.
logger = get_logger("main")


def error(msg):
    """
    The default error handler
    """
    global logger
    logger.error(f"ERROR: {msg}")
    sys.exit(1)
