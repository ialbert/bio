#
# Turn off broken pipe errors that may appear when piping into unix tools (head/tail etc)
#
# https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python/30091579#30091579
#
import sys
# Global package version
VERSION = "1.6.0"

try:
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)
except ImportError as exc:
    pass

# The minimal python version.
MIN_PY_VER = (3, 7)
if sys.version_info[:2] < MIN_PY_VER:
    msg1 = "# ERROR `bio` requires Python %i.%i or later.\n" % MIN_PY_VER
    msg2 = "# Your version appears to be %d.%d\n" % sys.version_info[:2]
    sys.stderr.write(msg1)
    sys.stderr.write(msg2)
    sys.exit(1)

