#
# Turn off broken pipe errors that may appear when piping into unix tools (head/tail etc)
#
# https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python/30091579#30091579
#
import sys

# Global package version
VERSION = "0.8.5"

try:
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)
except ImportError as exc:
    pass


