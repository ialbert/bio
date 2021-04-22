import sys, os

from biorun import utils

# Global alias file
ALIAS = dict()

alias_fname = os.path.join(utils.DATADIR, "alias.txt")

NOALIAS = '--noalias'

if NOALIAS in sys.argv:
    sys.argv.remove(NOALIAS)
elif  os.path.isfile(alias_fname):
    stream = open(alias_fname)
    stream = map(lambda x: x.strip().split(), stream)
    stream = filter(lambda x: len(x) > 1, stream)
    pairs = map(lambda x: x[:2], stream)
    ALIAS = dict(pairs)
