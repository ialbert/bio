"""
The main job runner. Register functions here.
"""
import importlib
import sys

import biorun.libs.placlib as plac
from biorun import utils, alias

# Module level logger
logger = utils.logger

#
# Valid subcommands:
#
# name = (module.function, auto_help, command_help)
#
SUB_COMMANDS = dict(
    fetch=("biorun.fetch.run", True, "downloads GenBank data from NCBI"),
    convert=("biorun.convert.run", True, "converts GenBank to FASTA or GFF"),
    meta=("biorun.meta.run", False, "downloads metadata by taxonomy ID"),
    taxon=("biorun.taxdb.run", False, "visualize taxonomies"),

    # define=("biorun.models.ontology.run", False, "explains biological terms"),
    # runinfo=("biorun.runinfo.run", True, "prints sequencing run information"),

)

# Generates the nicely indented help for each subcommand
block = [f"   bio {key:7} : {value[2]}" for (key, value) in SUB_COMMANDS.items()]

# Join help into a section.
block = "\n".join(block)

# Fill help section into the usage.
USAGE = f"""
bio: making bioinformatics fun again

Valid commands:

{block}

Examples:

    bio fetch NC_045512 MN996532 > genomes.gb 
    bio convert genomes.gb  --fasta
    bio convert genomes.gb  --type CDS --gff
    bio taxon 2697049 --lineage  
"""


def proof_reader(value):
    """
    Allows more error tolerant parameter input. Remaps -start to --start, --F to -F
    """

    # Longform parameter
    twodash = value.startswith("--")

    # Shortform parameter.
    onedash = value.startswith("-") and not twodash

    # Parameter is single letter
    oneletter = len(value.strip("-")) == 1

    try:
        # Can the value be converted to a number
        float(value)
        isnum = True
    except ValueError as exc:
        isnum = False

    # Single letter but two dashes. Drop a leading dash.
    if oneletter and twodash and not isnum:
        value = value[1:]

    # One dash but more more than one letter. Add a dash.
    if onedash and not oneletter and not isnum:
        value = f"-{value}"

    return value


def interrupt(func):
    """
    Intercept keyboard interrupts.
    """

    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except KeyboardInterrupt:
            sys.exit(0)

    return wrapper


@interrupt
def router():
    """
    Route the tasks based on subcommands parameters.
    """

    # Print usage when no parameters are passed.
    if len(sys.argv) == 1:
        print(USAGE)
        sys.exit(1)

    # Lowercase the subcommand.
    sys.argv[1] = sys.argv[1].lower()

    # Check the subcommand.
    cmd = sys.argv[1]

    # Maintain compatibility with a prior use case (convert == view).
    cmd = "convert" if cmd == "view" else cmd

    # Raise an error is not a valid subcommand.
    if cmd not in SUB_COMMANDS:
        print(USAGE, file=sys.stderr)
        logger.error(f"invalid command: {cmd}")
        sys.exit(-1)

    # Remove the command from the list.
    sys.argv.remove(cmd)

    # Allow multiple forms of parameters to be used.
    sys.argv = list(map(proof_reader, sys.argv))

    # Delegate to the imported method
    modfunc, flag, help = SUB_COMMANDS[cmd]

    # Add the help flag if no other information is present beyond command.
    if flag and len(sys.argv) == 1:
        sys.argv.append("-h")

    # Format: module.function
    mod_name, func_name = modfunc.rsplit(".", maxsplit=1)

    # Dynamic imports to allow other functionality
    # to work even when required for certain subcommands may be missing.
    mod = importlib.import_module(mod_name)

    # Get the function of the module
    func = getattr(mod, func_name)

    # Execute the function with plac.
    plac.call(func)


if __name__ == '__main__':
    router()
