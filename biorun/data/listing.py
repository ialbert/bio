"""
Lists the files in the database
"""
import plac, os, glob, gzip
from biorun import utils
from . import parse

# Get the logger information.
logger = utils.logger


def print_file_list():
    """
    Returns a list of the files in the data directory
    """
    pattern = os.path.join(os.path.join(utils.DATADIR, '*.gb.gz'))
    matched = glob.glob(pattern)
    for path in matched:
        fsize = utils.sizeof_fmt(os.path.getsize(path))
        base, fname = os.path.split(path)
        fname = fname.rsplit(".", maxsplit=2)[0]
        collect = parse.collect_metadata(path)
        title = collect.get("definition", "")
        project = collect.get("bioproject", ".       ")
        sample = collect.get("sample", ".       ")
        print (f"{fsize}\t{fname}\t{project}\t{sample}\t{title}")

@plac.flg('verbose', "verbose mode, progress messages printed")
def run(verbose=False):

    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    # Produce the file listing
    print_file_list()


