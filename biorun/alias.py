import sys, os

from biorun import utils

# Global alias file
ALIAS = dict()

def parse(fname):
    stream = open(fname)
    stream = map(lambda x: x.strip(), stream)
    stream = filter(lambda x: not x.startswith("#"), stream)
    stream = map(lambda x: x.strip().split(), stream)
    stream = filter(lambda x: len(x) > 1, stream)
    pairs = map(lambda x: x[:2], stream)
    data = dict(pairs)
    return data

# Potential alias files.
fname1 = "alias.txt"
fname2 = os.path.join(utils.DATADIR, fname1)

fnames = [fname1, fname2]
fnames = list(filter(os.path.isfile, fnames))

NOALIAS = '--noalias'

if NOALIAS in sys.argv:
    sys.argv.remove(NOALIAS)
elif fnames:
    ALIAS = parse(fnames[0])

# Print the alias upon executing the module directly.
if __name__ == '__main__':

    for key, value in ALIAS.items():
        print (f"{key}\t{value}")

