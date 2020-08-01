import os
import sys
from . import VERSION
import plac

# Commands names.
CONVERT, ALIGN = "convert", "align"

from .align import pairwise

# Enabled commands.
COMMANDS = {
    CONVERT: None,
    ALIGN: pairwise.run,
}

# Help printed when no command is passed.
USAGE = f"""
Bioinformatics utilities: {VERSION}

Usage: bio COMMAND 

Data commands:

    {CONVERT:10s} - convert biological data to other formats
    
Operations:

    {ALIGN:10s} - align sequences with different algorithms
  
Get more help on each command with:

    bio COMMAND -h 
"""

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
