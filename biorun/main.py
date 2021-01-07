"""
The main job runner. Register functions here.
"""
import sys, importlib
import biorun.libs.placlib as plac

from biorun import utils, const

# Module level logger
logger = utils.logger


def proofreader(value):
    """
    Bridges the gap between short and longforms. Allows the use of both by remapping to canonical forms.
    Remaps -start to --start, --F to -F
    """

    # Longform parameter
    long = value.startswith("--")

    # Shortform parameter.
    short = not long and value.startswith("-")

    # Is it one character.
    onechar = len(value.strip("-")) == 1

    try:
        # Can the value be converted to a number
        float(value)
        isnum = True
    except Exception as exc:
        isnum = False

    # Single character but long form. Drop a leading dash.
    if onechar and long and not isnum:
        value = value[1:]

    # Short form but more than one char. Add a dash.
    if short and not onechar and not isnum:
        value = f"-{value}"

    return value


def interrupt(func):
    """
    Keeps from raising tracebacks on keyboard interrupts.
    """
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except KeyboardInterrupt:
            sys.exit(0)
    return wrapper

@interrupt
def router(arglist=[]):
    """
    Routes the tasks based on incoming parameters.
    """

    # Allow multiple forms of parameters to be used.
    arglist = sys.argv = list(map(proofreader, sys.argv))

    # Delayed imports to allow other functionality to work even when some required libraries may be missing.

    # Check the presence of subcommands
    for cmd, mod in const.SUB_COMMANDS:

        # Found subcommand in arguments.
        if cmd in arglist:

            # Add the help flag if no other information is present beyond subcommand.
            if len(arglist) == 2:
                arglist.append("-h")

            # Import the module.
            lib = importlib.import_module(mod)

            # Execute the module.
            plac.call(lib.run)

            # Only one module may be called
            return

    # Default action, no subcommand was passed.

    from biorun import fetch

    # Add the help flag if no other information is present.
    if len(sys.argv) == 1:
        sys.argv.append("-h")

    # Delegate parameter parsing to converter.
    plac.call(fetch.run)


if __name__ == '__main__':
    router()
