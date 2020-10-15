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


def smartname(text):
    """
    Splits an accession number by colon into acc:name
    """
    pass

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

def first(value):
    """
   Returns the first element of a list.

    >>> first(["A", "B"])
    'A'
    """

    return value[0] if isinstance(value, list) else value

def flatten(value, sep="|"):
    """
    Flattens values that may be lists.

    >>> flatten(["A", "B"])
    'A|B'
    """

    return sep.join(map(str, value)) if isinstance(value, list) else value

def guess_type(path):
    """
    Attempts to guess a file type from an extension.
    """
    name, ext = os.path.splitext(path)
    ext = ext.lower()
    ftype = TYPE_BY_EXTENSION.get(ext, "")
    return ftype


def get_template(fname, dirname=__TMPL_DIR):
    """
    Loads and returns the content of a file.
    """
    path = os.path.join(dirname, fname)
    text = open(path).read()
    return text


def render_text(text, context={}):
    """
    Renders a template with a context.
    """
    tmpl = jinja2.Template(text, trim_blocks=True, lstrip_blocks=True, autoescape=False,
                           undefined=jinja2.StrictUndefined)
    text = tmpl.render(context)
    return text


def render_file(fname, context, dirname=__TMPL_DIR):
    """
    Renders a template from a file.
    """
    text = get_template(fname=fname, dirname=dirname)
    result = render_text(text, context)
    return result


def resolve_fname(acc, format='json'):
    """
    Resolve a file name given an accession number.
    """
    ext = format.lower()
    fname = f"{acc}.{ext}.gz"
    fname = os.path.join(DATADIR, fname)
    return fname


def sizeof_fmt(num, suffix=''):
    for unit in ['', 'K', 'M', 'G']:
        if abs(num) < 1024.0:
            return "%.0f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, '??', suffix)

@time_it
def save_stream(stream, fname, trigger=50000):
    """
    Write a input 'stream' as the fname filename.
    """
    # Save a stream into file and print progess.
    # Use a temporary file in case the process fails.
    # We don't want a failed cache file.
    tmp = tempfile.NamedTemporaryFile(mode="w+t")

    sequence = count(1)
    logger.info(f"streaming data to file")
    for index, line in zip(sequence, stream):
        if (index % trigger) == 0:
            logger.info(f"downloaded {index:,d} lines")
        tmp.write(line)

    logger.info(f"total of {index:,d} lines")

    # Not sure if this is needed. Can't hurt.
    tmp.flush()

    # Rewind temporary file to beginning
    tmp.seek(0)

    # Copy over the content from the temporary file to the final gzipped file destination.
    out = gzip.open(fname, 'wt')
    for line in tmp:
        out.write(line)
    out.close()
    tmp.close()

    # User friendly file size.
    fsize = sizeof_fmt(os.path.getsize(fname))
    logger.info(f"saved {fsize} data to {fname}")
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


# Initialize the logger.
logger = get_logger("main")


def error(msg):
    """
    The default error handler
    """
    global logger
    logger.error(f"ERROR: {msg}")
    sys.exit(1)
