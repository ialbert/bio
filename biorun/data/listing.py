"""
Lists the files in the database
"""
import os, glob, gzip, re, sys
from biorun import utils

# Get the logger information.
logger = utils.logger

def print_data_list():
    """
    Returns a list of the files in the data directory
    """
    pattern = os.path.join(os.path.join(utils.DATADIR, '*.json.gz'))
    matched = glob.glob(pattern)

    # Extract the definition from the JSON without parsing it.
    patt = re.compile(r'(definition\":\s*)(?P<value>\".+?\")')
    collect = []
    for path in matched:
        fsize = utils.sizeof_fmt(os.path.getsize(path))
        base, fname = os.path.split(path)
        fname = fname.rsplit(".", maxsplit=2)[0]

        # Parse the first N lines
        stream = gzip.open(path, 'rt') if path.endswith('gz') else open(path, 'rt')
        text = stream.read(1000)
        match = patt.search(text)

        title = match.group("value") if match else ''
        title = title.strip('", ')

        # Trim the title
        stitle = title[:100]
        stitle = stitle + "..." if len(title) != len(stitle) else stitle

        collect.append((str(fsize), f"{fname:10s}", stitle))

    collect = sorted(collect, key=lambda x: x[2])
    for row in collect:
        line = "\t".join(row)
        print(line)


