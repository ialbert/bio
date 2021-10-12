import codecs, csv, sys
from itertools import islice
from biorun import utils
from biorun.libs import placlib as plac

logger = utils.logger

def decode(text):
    """
    Recognize string encodings: \t etc
    """
    return codecs.decode(text, 'unicode_escape')


def print_metadata(terms, limit=None, header=False):

    if header:
        hdr = "accession species host date location isolate species_name".split()
        print (",".join(hdr))

    for term in terms:
        lines = get_metadata(term, limit=limit)
        tmp = next(lines)
        for line in lines:
            print(line)

# https://api.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/rest-api/

def get_metadata(taxid, limit=None, complete=True):
    """
    Returns all accessions
    """
    import requests

    url = f"https://api.ncbi.nlm.nih.gov/datasets/v1alpha/virus/taxon/{taxid}/genome/table"

    complete_only = "true" if complete else "false"
    params = {
        'format': 'csv',
        'refseq_only': "false",
        'complete_only': complete_only,
        'table_fields': [
            'nucleotide_accession',
            'species_tax_id',
            'host_tax_id',
            'collection_date',
            'geo_location',
            'isolate_name',
            'species_name',
        ]
    }

    conn = requests.get(url, stream=True, params=params)
    lines = conn.iter_lines()

    lines = islice(lines, limit)

    if conn.status_code != 200:
        msg = f"HTTP status code: {conn.status_code}"
        utils.error(msg)

    lines = map(decode, lines)

    return lines



@plac.flg('header', "print header", abbrev='H')
@plac.opt('limit', "download limit", abbrev='L')
def run(header=False, limit=None, *terms):
    limit = int(limit) if limit else None
    print_metadata(terms, limit=limit, header=header)


if __name__ == '__main__':
    # Bony fish: 117565

    plac.call(run)
