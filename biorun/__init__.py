import os

#
# Turn off broken pipe errors that may appear when piping into unix tools (head/tail etc)
#
# https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python/30091579#30091579
#
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

__CURRENT_DIR = os.path.dirname(__file__)

VERSION = "0.0.2"

