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


def match_srr(text):
    """
    Pattern for SRR numbers.
    """
    return bool(SRR.search(text))


# Documentation at https://www.ebi.ac.uk/ena/portal/api
def get_ena_fields(db='ena'):
    """
    Returns all valid ENA fields for a database
    """
    params = dict(dataPortal=db, result='read_run')
    stream = get_request(url=ENA_FIELDS, params=params, sep="\t")
    fields = [r['columnId'] for r in stream]
    fields.sort()
    return fields


def get_ncbi(text, db="protein", format="json"):
    drops = "statistics properties oslt".split()
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    params = dict(db=db, format=format, id=text)

    stream = get_request(url, params=params, bulk=True)


    data = json.loads(stream)

    data = data.get("result", {})

    collect = []
    for key in data:
        if key == 'uids':
            continue

        entry = data[key]
        for drop in drops:
            if drop in entry:
                del entry[drop]

        res = {}
        for key in sorted(entry.keys()):
            res[key] = entry[key]
        collect.append(res)

    if not collect:
        collect = [ dict(error="invalid genbank id", db=db, value=f"{text}")]
    return collect, None


def get_srr(text, all=False, sep=None):
    """
    Returns a list of SRR data.
    """
    url = ENA_REPORT
    if all:
        fields = get_ena_fields()
    else:
        fields = [
            'run_accession',
            "sample_accession",
            'first_public',
            'country',
            'sample_alias',
            'fastq_bytes',
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
    )

    stream = get_request(url, params=params, sep=sep)

    return stream, None


def match_bioproject(text):
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
def match_genbank_nucleotide(text):
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


def match_genbank_protein(text):
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


def dictreader(stream, sep=None):
    """
    Function to wrap a stream into a DictReader.
    """
    return csv.DictReader(stream, delimiter=sep)


def get_request(url, params={}, sep=None, bulk=False):
    try:

        # print (url, params, file=sys.stderr)

        r = requests.get(url, params=params)

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        if bulk:
            return r.text
        else:
            stream = r.iter_lines(decode_unicode=True)
            stream = dictreader(stream, sep=sep) if sep else stream
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

    if len(hits) < total:
        warn = f'#  showing {len(hits)} out of {total} results.'
    else:
        warn = None

    return hits, warn


def dispatch(word, all=False, fields='', limit=5, species='', scopes=''):
    if match_srr(word) or match_bioproject(word):
        values, warn = get_srr(word, all=all, sep="\t")

    elif match_genbank_nucleotide(word):
        values, warn = get_ncbi(word, db="nuccore")

    elif match_genbank_protein(word):
        values, warn = get_ncbi(word, db="protein")

    else:
        fields = ",".join(['symbol', 'name', 'taxid', fields])
        values, warn = search_mygene(word, fields=fields, limit=limit, species=species, scopes=scopes)

    return values, warn


@plac.flg('csv_', "produce comma separated output")
@plac.flg('tab', "produce tab separated output")
@plac.flg('all', "get all possible fields")
@plac.flg('header', "show headers", abbrev="H")
@plac.opt('limit', "download limit", abbrev='l')
@plac.opt('fields', "fields", abbrev='f')
@plac.opt('species', "species", abbrev='s')
@plac.opt('scopes', "scopes", abbrev='S')
@plac.pos('query', "query terms")
def run(all=False, csv_=False, tab=False, header=False, species='', scopes='symbol', limit=5, fields='', *words):


    sep = None

    sep = "," if csv_ else sep

    sep = "\t" if tab else sep

    collect = []
    warns = []
    for word in words:
        values, warn = dispatch(word, all=all, limit=limit, fields=fields, species=species, scopes=scopes)

        collect.extend(values)
        warns.append(warn)

    if sep:
        fields = collect[0].keys()
        wrt = csv.writer(sys.stdout, delimiter=sep)
        stream = [ x.values() for x in collect]
        keys = [ x.keys() for x in collect]
        if header:
            wrt.writerow(keys[0])
        wrt.writerows(stream)
    else:
        text = json.dumps(collect, indent=4)
        print(text)

    # Show collected warnings at the end where it is not so easy to miss.
    warns = filter(None, warns)
    for warn in warns:
        print(warn, file=sys.stderr)


# SRR5260547
if __name__ == '__main__':
    plac.call(run)
