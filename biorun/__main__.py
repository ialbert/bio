"""
The main job runner. Register functions here.
"""
import os
import sys
from biorun import VERSION
import plac
from biorun import utils

# Commands names.
CONVERT, ALIGN, FETCH = "convert", "align", "fetch"

from biorun.align import pairwise
from biorun.data import fetch

# Enabled commands.
COMMANDS = {
    CONVERT: None,
    FETCH: fetch.run,
    ALIGN: pairwise.run,
}

# Context for the USAGE help page.
context = dict(VERSION=VERSION,
               CONVERT=CONVERT, ALIGN=ALIGN, FETCH=FETCH
               )

# The default help page for the tool.
USAGE = utils.render_file("usage.txt", context=context)

@plac.annotations(
    cmd="command"
)
def run(*cmd):
    """
    Main command runner.
    """

    # Commands are case insensitive.
    target = cmd[0].lower() if cmd else None

    # Print usage if no command is seen.
    if target is None or target in ("-h", "--help"):
        print(f"{USAGE}", file=sys.stderr)
        sys.exit(1)

    # Handles invalid command.
    if target not in COMMANDS:
        print(f"{USAGE}", file=sys.stderr)
        print(f"*** Invalid command: {target}\n", file=sys.stderr)
        sys.exit(127)

    try:
        # Call the targetet function.
        func = COMMANDS[target]
        plac.call(func, cmd[1:])
    except KeyboardInterrupt:
        # Breakout from interrupts without a traceback.
        sys.exit(0)

run.add_help = False

def main():
    plac.call(run)

if __name__ == '__main__':
    main()
