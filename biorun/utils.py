"""
Utilites funcions.
"""
import os
import sys
from itertools import count, islice
import time
import jinja2

from biorun import DUMP_DIR

# The path to the current file.
__CURR_DIR = os.path.dirname(__file__)

# The default path to templates.
__TMPL_DIR = os.path.join(__CURR_DIR, "templates")

OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'

GENBANK, FASTA, GFF, BED, SAM, BAM = "genbank", "fasta", "gff", "bed", "sam", "bam"

TYPE_BY_EXTENSION = {
    "gb": GENBANK, "gbk": GENBANK, "genbank": GENBANK,
    "fa": FASTA, "fasta": FASTA,
    "bed": BED,
    "gff": GFF,
    "sam": SAM,
    "bam": BAM,
}

def guess_type(path):
    """
    Attempts to guess a file type from an extension.
    """
    name, ext = os.path.splitext(path)
    ext = ext.lower()
    ftype = TYPE_BY_EXTENSION.get(ext, "")
    return ftype


def print_message(styles=[OKBLUE], msg='', verb=0):
    styles = ''.join(styles)
    if verb >= 1:
        print(f"{styles}{msg}{ENDC}")


def get_template(fname, dirname=__TMPL_DIR):
    """
    Loads and returns the content of a file.
    """
    path = os.path.join(dirname, fname)
    text = open(path).read()
    return text


def timer(func):
    """
    Decorator used to time functions.
    """
    def __wrapper__(*args, **kwargs):

        t1 = time.time()
        res = func(*args, **kwargs)
        delta = time.time() - t1
        print(f"{func.__name__} : {delta} seconds.")
        return res
    return __wrapper__


def timer_func(verb):
    """
    Prints progress on inserting elements.
    """

    last = time.time()

    def elapsed(msg):
        nonlocal last
        now = time.time()
        sec = round(now - last, 1)
        last = now
        print_message(styles=[OKGREEN], msg=f"{msg} in {sec} seconds.", verb=verb)

    def progress(index, step=5000, msg=""):
        nonlocal last
        if index % step == 0:
            elapsed(f"... {index} {msg}")

    return elapsed, progress


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


def resolve_fname(acc, directory=None, ext='gbk'):
    """
    Resolve a file name given an accession number.
    """
    suffix = f"{acc}.{ext}"
    directory = directory or DUMP_DIR
    fname = os.path.abspath(os.path.join(directory, suffix))
    return fname


def download(stream, outname, buffer=1024, verb=0, formatter=lambda x: x):
    """
    Write a input 'stream' into the output filename.
    Overwrite existing file given a flag.
    """

    # Ensure directory exists.
    outdir = os.path.dirname(outname)
    os.makedirs(outdir, exist_ok=True)
    stream = islice(zip(count(1), stream), None)
    elapsed, progress = timer_func(verb=verb)

    # Write 'stream' into output and print progess
    with open(outname, 'w', buffering=buffer) as output_stream:
        for index, line in stream:
            progress(index, msg="lines", step=500)
            output_stream.write(formatter(line))

    # Show file size in Mb
    fsize = os.path.getsize(outname) / 1000 / 1000
    elapsed("Wrote {0:1f} MB to {1}".format(fsize, outname))
    return

