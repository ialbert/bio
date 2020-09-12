"""
The main job runner. Register functions here.
"""
import os, time
import sys
import plac
from importlib import import_module
from biorun import VERSION
from biorun import utils

# Commands names.
CONVERT, ALIGN, FETCH = "convert", "align", "fetch"

# Enabled commands.
COMMANDS = {
    CONVERT: 'biorun.data.convert',
    FETCH: 'biorun.data.fetch',
    ALIGN: 'biorun.align.pairwise'
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
        # Import and call the target function.
        package = COMMANDS[target]
        module = import_module(name=package)
        # Target the run function in each module.
        func = getattr(module, 'run')
        rest = cmd[1:]
        rest = rest or [ "-h"]
        plac.call(func, rest)
    except KeyboardInterrupt:
        # Terminate keyboard interrupts without a traceback.
        sys.exit(0)


run.add_help = False

def main():
    plac.call(run)

if __name__ == '__main__':
    main()
