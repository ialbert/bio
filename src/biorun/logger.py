import logging

def get_logger(name="bio", hnd=None, fmt=None, terminator='\n'):
    """
    Initializes a logger with a handler and formatter.
    """
    # Get the logger name.
    log = logging.getLogger(name)

    # Default logging level.
    log.setLevel(logging.WARNING)

    # The log handler.
    hnd = hnd or logging.StreamHandler()
    hnd.terminator = terminator

    # The logging formatter.
    fmt = fmt or logging.Formatter('# %(message)s')

    # Add formatter to handler
    hnd.setFormatter(fmt)

    # Add handler to logger
    log.addHandler(hnd)

    return log