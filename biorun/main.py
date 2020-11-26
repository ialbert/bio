"""
The main job runner. Register functions here.
"""
import sys
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
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except KeyboardInterrupt:
            sys.exit(0)
    return wrapper

@interrupt
def router():
    """
    Routes the tasks based on incoming parameters.
    """

    # Allow multiple forms of parameters to be used.
    sys.argv = list(map(proofreader, sys.argv))

    # Delayed imports to allow other functionality to work even when some required libraries may be missing.

    if const.ALIGN_COMMAND in sys.argv:
        # Run alignment related functionality..
        sys.argv.remove(const.ALIGN_COMMAND)
        from biorun.methods import align
        plac.call(align.run)

    elif const.TAXON_COMMAND in sys.argv:
        # Run taxonomy related functionality.
        from biorun.models import taxdb
        sys.argv.remove(const.TAXON_COMMAND)
        plac.call(taxdb.run)

    elif const.DBLINK_COMMAND in sys.argv:
        # Run SRA specific functionality.
        from biorun.models import dblink
        sys.argv.remove(const.DBLINK_COMMAND)
        plac.call(dblink.run)

    else:
        # Default action is to convert a file.
        from biorun import convert

        # Add the help flag if no other information is present.
        if len(sys.argv) == 1:
            sys.argv.append("-h")

        # Delegate parameter parsing to converter.
        plac.call(convert.run)


if __name__ == '__main__':
    router()
