"""
Partially parses GenBank files.

Used when only a certain small/tiny bit of information needs to be extracted that does not warrant fully
parsing the entire potentially huge file.
"""
import gzip, io
from biorun import models
from itertools import *

def open_stream(path):
    if path.endswith(".gz"):
        return gzip.open(path, 'rt')
    else:
        return open(path, 'rt')

def collect_metadata(path, limit=100):
    """
    Attempts to find simple tidbits from a genbank file
    """
    collect = dict()
    stream = open_stream(path)
    func = lambda x: not x.startswith("FEATURES")
    stream = takewhile(func, stream)
    for line in stream:
        if line.startswith("DEFINITION"):
            value = line.split(maxsplit=1)[1].strip()
            collect["definition"] = value.strip()
        if "BioProject:" in line:
            value = line.split()[-1].strip()
            collect["bioproject"] = value
        if "BioSample:" in line:
            value = line.split()[-1].strip()
            collect["sample"] = value

    return collect