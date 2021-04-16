import codecs
from itertools import islice
from biorun import utils
from biorun.libs import placlib as plac

logger = utils.logger

def decode(text):
    """
    Recognize string encodings: \t etc
    """
    return codecs.decode(text, 'unicode_escape')


def print_metadata(terms, limit=None):
    def formatter(row):
        print("\t".join(row))

    for term in terms:
        lines = get_metadata(term, limit=limit)
        lines = filter(lambda x: x.split(), lines)
        old_header = next(lines)
        new_header = "accession species host date location isolate species_name".split()

        print("\t".join(new_header))
        for line in lines:
            print(line)


def get_metadata(taxid, limit=None, complete=True):
    """
    Returns all accessions
    """
    import requests

    url = f"https://api.ncbi.nlm.nih.gov/datasets/v1alpha/virus/taxon/{taxid}/genome/table"

    complete_only = "true" if complete else "false"
    params = {
        'format': 'tsv',
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



@plac.flg('taxon', "downloads metadata for the taxon", abbrev='m')
@plac.flg('viral', "downloads viral metadata", abbrev='v')
@plac.opt('limit', "download limit", abbrev='L')
def run(taxon=False, viral=False, limit=None, *terms):
    limit = int(limit) if limit else None
    print_metadata(terms, limit=limit)


if __name__ == '__main__':
    # Bony fish: 117565

    plac.call(run)
