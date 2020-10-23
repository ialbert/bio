"""
Lists the files in the database
"""
import plac, os, glob, gzip, re, sys
from biorun import utils

# Get the logger information.
logger = utils.logger

def print_file_list():
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

@plac.opt('alias', "create alias to accession ")
@plac.flg('delete', "delete the data for accession")
@plac.flg('verbose', "verbose mode, progress messages printed")
def run(alias='', delete=False, verbose=False, *accs):
    # Set the verbosity of the process.
    utils.set_verbosity(logger, level=int(verbose))

    for acc in accs:

        # Make aliases
        if alias:

            json_name1 = utils.resolve_fname(acc=acc, format="json")
            json_name2 = utils.resolve_fname(acc=alias, format="json")

            names = [
                (json_name1, json_name2),
            ]
            for src, dest in names:
                utils.symlink(src, dest)
                logger.info(f"link {dest}")

        elif delete:

                json_name = utils.resolve_fname(acc=acc, format="json")
                if os.path.isfile(json_name):
                    os.remove(json_name)


    # Produce the file listing
    print_file_list()
