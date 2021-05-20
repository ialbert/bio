import os
import sys

from biorun import utils

# Global alias file
ALIAS = dict()


def parse(fname):
    """
    Parses an alias file, returns a dictionary mapping: Accession -> Name
    """
    stream = open(fname)
    stream = map(lambda x: x.strip(), stream)
    stream = filter(lambda x: not x.startswith("#"), stream)
    stream = map(lambda x: x.strip().split(), stream)
    stream = filter(lambda x: len(x) > 1, stream)
    pairs = map(lambda x: x[:2], stream)
    data = dict(pairs)
    return data


def safe_parse(fname):
    """
    Don't raise errors on invalid alias files.
    """
    try:
        return parse(fname)
    except Exception as exc:
        return {}


# Alias files can be local or global.
alias_local = "alias.txt"
alias_global = os.path.join(utils.DATADIR, alias_local)

# Find the first existing alias file.
fnames = [alias_local, alias_global]
fnames = list(filter(os.path.isfile, fnames))

# This flag turns off aliasing.
NOALIAS = '--noalias'

# Remove the noalias if passed.
if NOALIAS in sys.argv:
    sys.argv.remove(NOALIAS)
elif fnames:
    # Attempts to read the alias file.
    ALIAS = safe_parse(fnames[0])

# Print the alias upon executing the module directly.
if __name__ == '__main__':

    for key, value in ALIAS.items():
        print(f"{key}\t{value}")
