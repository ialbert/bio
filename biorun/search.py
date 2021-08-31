import json

import plac


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
def runme(query, fields, species='', scopes='', size=3):
    import mygene

    mg = mygene.MyGeneInfo()
    # res = mg.querymany(names, scopes='symbol', fields='symbol,name')

    data = mg.query(query, fields=fields, scopes=scopes, species=species, size=size)

    text = json.dumps(data, indent=4)

    print(text)


@plac.flg('header', "print header", abbrev='H')
@plac.opt('limit', "download limit", abbrev='l')
@plac.opt('fields', "fields", abbrev='f')
@plac.opt('species', "species", abbrev='s')
@plac.opt('scopes', "scopes", abbrev='S')
@plac.pos('query', "download limit")
def run(header=False, limit=5, species='', fields='', scopes='symbol', *query):
    limit = int(limit) if limit else None


    fields = ",".join(['symbol', 'name', 'taxid', fields])

    query = ' '.join(query)

    runme(query, fields=fields, species=species, scopes=scopes, size=limit)


if __name__ == '__main__':
    # Bony fish: 117565

    plac.call(run)
