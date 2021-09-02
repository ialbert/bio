import json

import plac

import mygene

#
# Available fields
#
# https://docs.mygene.info/en/latest/doc/query_service.html#available_fields
#
# Guide to jq
#
# https://earthly.dev/blog/jq-select/
#
# bio gene symbol:cdk2 -limit 1 -fields refseq -species 9823
#
# cat out.json | jq -r '.hits[].refseq.translation[]|[.protein, .rna] | @tsv'
#
from biorun import taxon

def execute(query, fields, species='', scopes='', size=3):

    client = mygene.MyGeneInfo()

    data = client.query(query, fields=fields, scopes=scopes, species=species, size=size)

    # Fill in taxonomy name to the database
    names, graph = taxon.get_data(strict=False)

    for hit in data.get('hits', []):
        hit['taxname'] = names.get(hit.get('taxid'), [''])[0]
    text = json.dumps(data, indent=4)

    print(text)


@plac.opt('limit', "download limit", abbrev='l')
@plac.opt('fields', "fields", abbrev='f')
@plac.opt('species', "species", abbrev='s')
@plac.opt('scopes', "scopes", abbrev='S')
@plac.pos('query', "download limit")
def run(limit=5, species='', fields='', scopes='symbol', *query):
    limit = int(limit) if limit else None

    fields = ",".join(['symbol', 'name', 'taxid', fields])

    query = ' '.join(query)

    execute(query, fields=fields, species=species, scopes=scopes, size=limit)


if __name__ == '__main__':
    # Bony fish: 117565

    plac.call(run)
