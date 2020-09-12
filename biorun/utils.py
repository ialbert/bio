"""
Utilites funcions.
"""
import sys, os, tempfile, shutil
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
        start = time.time()
        try:
            return func(*args, **kwargs)
        finally:
            end = time.time()
            diff = int(round((end-start), 1)) or 0.1
            logger.info(f"{func.__name__} execution time: {diff} seconds")
    return timer

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


def resolve_fname(acc, format='gb'):
    """
    Resolve a file name given an accession number.
    """
    ext = format.lower()
    ext = 'fa' if ext == 'fasta' else ext
    fname = f"{acc}.{ext}"
    fname = os.path.join(DATADIR, fname)
    return fname


def save_stream(stream, fname, trigger=100000, stdout=False):
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
            logger.info(f"processed {index:,d} lines")
        tmp.write(line)
        if stdout:
            print(line, end='')

    # Not sure if this is needed. Can't hurt.
    tmp.flush()

    # Rewind temporary file to beginning
    tmp.seek(0)

    # Copy over the content from the temporary file to the final destination
    res = open(fname, 'wt')
    for line in tmp:
        res.write(line)
    res.close()
    tmp.close()

    # Show file size in Mb
    fsize = os.path.getsize(fname) / 1024 / 1024
    logger.info(f"saved {fsize:0.1f} MB data to {fname}")
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
