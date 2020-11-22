"""
Attempts to produce the SRA and GEO based information for the record.

The information is not always properly documented in the genbank file.
"""
from biorun.libs import placlib as plac

@plac.pos("acc", "accessions")
@plac.flg('inter', "interactive (data from command line)")
@plac.flg('verbose', "verbose mode, prints more messages")

def run(inter=False, verbose=False, **acc):
    print ("SRA")
    pass
