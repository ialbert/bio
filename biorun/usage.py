#
# Help template for the usage
#

HELP = '''
Bioinformatics utilities: {VERSION}

Usage: bio COMMAND

Data commands

   {LIST:10s} - lists the available data in the cache
   {FETCH:10s} - downloads GenBank/EMBL data into cache
   {VIEW:10s} - converts/extracts information from GenBank/EMBL data

Action commands

   {ALIGN:10s} - aligns sequences with different algorithms

Get more help on each command with:

    bio COMMAND -h

Examples:

    bio view AF086833
'''