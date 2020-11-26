"""
Attempts to produce the SRA and GEO based information for the record.

The information is not always properly documented in the genbank file.
"""
from biorun.libs import placlib as plac
from biorun import storage
from Bio import Entrez

Entrez.email = "bio@example.org"

PROJECT, SAMPLE = "BioProject:", "BioSample:"

def search(term):
    handle = Entrez.esearch(db='sra', term=term)
    record = Entrez.read(handle)
    handle.close()
    print (record)

def parse_data(rec):
    """
    Attempts to extract information from a record.
    """
    store = dict()
    dblinks = rec.get("dblink", [])
    for link in dblinks:
        key, value = link.split(":")
        store[key] = value

    return store

def process_storage(acc):

    collect = []

    for name in acc:
        data = storage.get_json(name) or []
        for rec in data:
            store = parse_data(rec)
            collect.append(store)

    return collect

def print_dict(rec):

    for key, value in rec.items():
        print (f"{key}\t{value}")

def get_runinfo(term):
    pass

@plac.pos("acc", "accessions")
@plac.flg('inter', "interactive (data from command line)")
@plac.flg('verbose', "verbose mode, prints more messages")
def run(inter=False, verbose=False, *acc):

    terms = acc if inter else process_storage(acc)

    terms = list(filter(None, terms))

    for rec in terms:
        print_dict(rec)

    #for term in terms:
    #    search(term)

    pass
