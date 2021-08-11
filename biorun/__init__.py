#
# Turn off broken pipe errors that may appear when piping into unix tools (head/tail etc)
#
# https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python/30091579#30091579
#
import sys

# Global package version
VERSION = "0.6.2"

try:
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)
except ImportError as exc:
    pass

# Import here to detect missing installation early on.
try:
    from Bio import Entrez
except ImportError as exc:
    print(f"### Error: {exc}", file=sys.stderr)
    print(f"### This program requires biopython", file=sys.stderr)
    print(f"### Install: conda install -y biopython>=1.78", file=sys.stderr)
    sys.exit(-1)


