import sys, plac
from biorun import utils

@plac.pos("seq", "sequences")
@plac.opt("sea", "database", choices=["nuccore", "protein"])
@plac.opt("format_", "return format", choices=["gbwithparts", "fasta", "gb"])
def run(db="nuccore", format_="gbwithparts", alias='', *acc):

    ids = []
    for num in acc:
        ids.extend(num.split(","))

    if not sys.stdin.isatty():
        ids.extend(utils.read_lines(sys.stdin))

if __name__ == '__main__':
    # id = "AY851612",

    run()
