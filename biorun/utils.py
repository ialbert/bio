"""
Utilites funcions.
"""
import os
import jinja2

# The path to the current file.
__CURR_DIR = os.path.dirname(__file__)

# The default path to templates.
__TMPL_DIR = os.path.join(__CURR_DIR, "templates")


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

