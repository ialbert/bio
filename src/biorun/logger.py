import logging
import sys

DEBUG = logging.DEBUG
WARNING = logging.WARNING
INFO = logging.INFO
ERROR = logging.ERROR

def get_logger(name="bio", hnd=None, fmt=None, terminator='\n'):
    """
    Initializes a logger with a handler and formatter.
    """
    # Get the logger name.
    log = logging.getLogger(name)

    # Default logging level.
    log.setLevel(logging.WARNING)

    # The log handler.
    hnd = hnd or logging.StreamHandler(sys.stderr)
    hnd.terminator = terminator

    # The logging formatter.
    fmt = fmt or logging.Formatter('# %(levelname)s\t%(module)s.%(funcName)s\t%(message)s')

    # Add formatter to handler
    hnd.setFormatter(fmt)

    # Add handler to logger
    log.addHandler(hnd)

    return log

# Default logger.
default = get_logger("main")