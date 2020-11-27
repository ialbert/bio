"""
Attempts to produce the SRA and GEO based information for the record.

The information is not always properly documented in the genbank file.
"""
from biorun.libs import placlib as plac
from biorun import storage, utils
from Bio import Entrez
from biorun import entrez
import json, sys, csv

Entrez.email = "bio@example.org"

PROJECT, SAMPLE = "BioProject", "BioSample"

# The common logger
logger = utils.logger

def pprint(data):
    print(json.dumps(data, indent=4))


def search(term, db='sra', tabular=False):

    env = entrez.esearch(db=db, term=term, usehistory="y")

    data = entrez.efetch(db=db, env=env, retmax=2, rettype="runinfo")

    elems = data['SraRunInfo']["Row"]

    if tabular and elems:
        fieldnames = elems[0].keys()
        writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in elems:
            #print (row)
            writer.writerow(row)

    else:
        pprint(elems)

    return data



def parse_dblinks(rec):
    """
    Attempts to extract information from a record.
    """
    store = {PROJECT: '', SAMPLE: ''}
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
            store = parse_dblinks(rec)
            collect.append((name, store))

    return collect


def print_links(pair):
    """
    Prints the content of the dblinks.
    """
    acc, data = pair
    for key, value in data.items():
        print(f"{acc}\t{key}\t{value}")


def get_runinfo(term):
    pass


@plac.pos("acc", "accessions")
@plac.flg('inter', "interactive (data from command line)")
@plac.flg('project', "project run information")
@plac.flg('sample', "sample information")
@plac.flg('table', "tabular output")
@plac.flg('verbose', "verbose mode, prints more messages")
def run(project=False, sample=False, table=False, inter=False, verbose=False, *acc):

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if inter:
        # The query terms will be the same as the input
        collect = [ (t, {PROJECT:t, SAMPLE:t}) for t in acc]
    else:
        # Parse the query terms from the data
        collect = process_storage(acc)

    for row in collect:
        name, metadata = row
        if project:
            term = metadata[PROJECT]
            search(term, tabular=table)
        elif sample:
            term = metadata[SAMPLE]
            search(term, tabular=table)
        else:
            print_links(row)

