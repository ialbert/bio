import csv
import json
import re
import sys

import requests

from biorun import utils
from biorun.libs import placlib as plac

# SRR numbers: SRR5260547
SRR = re.compile(r'(ERR|SRR)\d+')

# Bioproject numbers: PRJNA374918
PRJ = re.compile(r'PRJ([A-Z])+\d+')

# Genbank accessions: NC_045512
GBK = re.compile(r'(?P<letters>[a-zA-Z]+)(?P<under>_?)(?P<digits>\d+)(\.(?P<version>\d+))?')

# ENA API points
ENA_API = "https://www.ebi.ac.uk/ena/portal/api"
ENA_FIELDS = f"{ENA_API}/returnFields"
ENA_REPORT = f"{ENA_API}/filereport"


def is_srr(text):
    """
    Pattern for SRR numbers.
    """
    return bool(SRR.search(text))


# Documentation at https://www.ebi.ac.uk/ena/portal/api

def get_srr(text, all=False, sep=False, limit=100):
    """
    Performs an SRR search
    """
    url = ENA_REPORT
    if all:
        params = dict(dataPortal='ena', result='read_run')
        stream = get_request(ENA_FIELDS, params=params, sep="\t")
        fields = [r['columnId'] for r in stream]
        fields.sort()
    else:
        fields = [
            'run_accession',
            "sample_accession",
            'first_public',
            'country',
            'sample_alias',
            'read_count',
            'library_name',
            "library_strategy",
            "library_source",
            'library_layout',
            'instrument_platform', 'instrument_model',
            'study_title',
            'fastq_ftp',
        ]

    fields = ",".join(fields)

    params = dict(
        accession=text,
        fields=fields,
        result='read_run',
        #limit=limit,
    )
    stream = get_request(url, params=params, sep="\t")

    if sep:

        wrt = csv.DictWriter(sys.stdout, fieldnames=stream.fieldnames, delimiter=sep)
        wrt.writerows(stream)
    else:
        stream = list(stream)
        text = json.dumps(stream, indent=4)
        print(text)


def is_bioproject(text):
    """
    Pattern for project numbers.
    """
    return bool(PRJ.search(text))


def parse_genbank(text):
    """
    Attempts to parse text into a NCBI structure.
    """
    m = GBK.search(text)
    code = m.group("letters") if m else ''
    digits = m.group("digits") if m else ''
    refseq = m.group("under") if m else ''
    version = m.group("version") if m else ''
    return code, digits, refseq, version


#
# https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
#
def is_genbank_nucleotide(text):
    """
    Returns true if text matches NCBI nucleotides.
    """
    code, digits, refseq, version = parse_genbank(text)
    if refseq:
        cond = code in ["AC", "NC", "NG", "NT", "NW", "NZ", "NM", "XM", "XR"]
    else:
        num1, num2 = len(code), len(digits)
        cond = (num1 == 1 and num2 == 5) or (num1 == 2 and num2 == 6) or (num1 == 3 and num2 == 8)
    return cond


def is_genbank_protein(text):
    """
    Returns true if text matches NCBI protein sequences
    """
    code, digits, refseq, version = parse_genbank(text)
    if refseq:
        cond = code in ["AP", "NP", "YP", "XP", "WP"]
    else:
        num1, num2 = len(code), len(digits)
        cond = (num1 == 3 and num2 == 5) or (num1 == 3 and num2 == 7)
    return cond


def get_request(url, params={}, sep=None, bulk=False):
    try:

        #print (url, params, file=sys.stderr)

        r = requests.get(url, params=params)

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        if bulk:
            return r.text
        else:
            stream = r.iter_lines(decode_unicode=True)
            if sep:
                stream = csv.DictReader(stream, delimiter=sep)
            return stream

    except Exception as exc:
        utils.error(f"Error for {url}, {params}: {exc}")



def search_mygene(query, fields, species='', scopes='', limit=5):
    import mygene
    from biorun import taxon

    client = mygene.MyGeneInfo()

    data = client.query(query, fields=fields, scopes=scopes, species=species, size=limit)

    total = data.get('total', 0)

    # Get rid of data we don't need
    hits = data.get('hits', [])

    # Fill in taxonomy name to the database
    names, graph = taxon.get_data(strict=False)

    # Fill the taxonmy name, get rid of fields we don't want.
    for hit in hits:
        del hit['_id']
        del hit['_score']
        hit['taxname'] = names.get(hit.get('taxid'), [''])[0]

    text = json.dumps(hits, indent=4)

    print(text)

    if len(hits) < total:
        print(f'#  showing {len(hits)} out of {total} results.', file=sys.stderr)


def dispatch(word, all=False, sep=None, fields='', limit=5, species='', scopes=''):

    if is_srr(word) or is_bioproject(word):
        get_srr(word, all=all, sep=sep)
    else:
        fields = ",".join(['symbol', 'name', 'taxid', fields])
        search_mygene(word, fields=fields, limit=limit, species=species, scopes=scopes)

@plac.flg('csv', "produce comma separated output")
@plac.flg('tab', "produce tab separated output")
@plac.flg('all', "get all possible fields")
@plac.opt('limit', "download limit", abbrev='l')
@plac.opt('fields', "fields", abbrev='f')
@plac.opt('species', "species", abbrev='s')
@plac.opt('scopes', "scopes", abbrev='S')
@plac.pos('query', "query terms")
def run(all=False, csv=False, tab=False, species='', scopes='symbol', limit=5, fields='', *words):

    sep = None

    sep = "," if csv else sep

    sep = "\t" if tab else sep

    for word in words:
        dispatch(word, all=all, sep=sep, limit=limit, fields=fields, species=species, scopes=scopes)


# SRR5260547
if __name__ == '__main__':
    plac.call(run)
