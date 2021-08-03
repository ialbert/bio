"""
Generate command lines to obtain data from NCBI.
"""
import os, sys, re
from itertools import *
from collections import namedtuple
from urllib.request import urlopen

import biorun.libs.placlib as plac
from biorun import utils


# Guess accession numbers that are proteins.
# https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/

# Ncbi accession pattern
ncbi_patt = r'(?P<letters>[a-zA-Z]+)(?P<under>_?)(?P<digits>\d+)(\.(?P<version>\d+))?'
ncbi_patt = re.compile(ncbi_patt)

NcbiAcc = namedtuple("Accesion", "value letters digits refseq version")


# Matching hrefs in pages.
href_patt = r'href=[\'"]?([^\'" >]+)'
href_patt = re.compile(href_patt)

def error(msg, stop=True):
    print (f"# *** {msg}", file=sys.stderr)
    if stop:
        sys.exit(1)

def detect_format(text):
    m = ncbi_patt.search(text)
    if not m:
        error(f"accession number format not recognized: {text}")

    letters, digits, under, version = m.group("letters"), m.group("digits"), m.group("under"), m.group("version")

    acc = NcbiAcc(value=text, letters=letters, digits=digits, refseq=bool(under), version=version)

    return acc

# Refseq accession numbers
# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/

def is_nucleotide(acc):
    if acc.refseq:
        return acc.letters in ["AC", "NC", "NG", "NT", "NW", "NZ", "NM", "XM", "XR"]

    num1, num2 = len(acc.letters), len(acc.digits)
    cond = (num1==1 and num2 ==5) or (num1==2 and num2==6) or (num1==3 and num2==8)

    return cond


def is_protein(acc):
    if acc.refseq:
        return acc.letters in [ "AP", "NP", "YP", "XP", "WP"]

    num1, num2 = len(acc.letters), len(acc.digits)
    cond = (num1 == 3 and num2 == 5) or (num1 == 3 and num2 == 7)

    return cond


def is_assembly(acc):
    return acc.letters == "GCF" or acc.letters == "GCA"

def make_csv(accs):
    return ",".join([acc.value for acc in accs])

def read_url(url):
    try:
        text = urlopen(url).read().decode("utf-8")
    except Exception as exc:
        msg = f"error reading: {url}"
        error(msg, stop=False)
        text = ''

    return text

def get_links(text):

    matches = href_patt.findall(text)
    for match in matches:
        if not match.startswith("/"):
            yield match
    return text


def write_command(ids, db, format_):
    #ids = ",".join(ids)

    accs = list(map(detect_format, ids))

    # Sanity check
    nucs = [ acc for acc in accs if is_nucleotide(acc) ]
    prot = [ acc for acc in accs if is_protein(acc)]
    asem = [ acc for acc in accs if is_assembly(acc)]

    n_vals, p_vals, a_vals = make_csv(nucs), make_csv(prot), make_csv(asem)

    if nucs and prot:
        error(f"nucl = {n_vals}", stop=False)
        error(f"prot = {p_vals}", stop=False)
        error("may not mix nucleotide and protein accessions")

    if asem and nucs:
        error(f"nucl = {n_vals}", stop=False)
        error(f"asem = {a_vals}", stop=False)
        error("may not mix assembly accessions with nucleotide accessions")

    if asem and prot:
        error(f"prot = {p_vals}", stop=False)
        error(f"asem = {a_vals}", stop=False)
        error("may not mix assembly accessions with protein accessions")

    if prot:
        db = "protein"

    if not asem:
        line = ",".join(ids)
        cmd = f"epost -db {db} -id {line} | efetch -format {format_}"
        print(cmd)
        return

    # Generate each report
    for elem in asem:

        value = elem.value
        print(f"# {value}")
        # Figure out if the assembly is present online or not.

        root = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
        p1, p2, p3, = value[4:7], value[7:10], value[10:13]
        url1 = f"{root}/{elem.letters}/{p1}/{p2}/{p3}/"

        print (f"# {url1}")
        text = read_url(url1)
        if not text:
            cmd = f"esearch -db assembly -query {value} | elink -target nuccore | efetch -format fasta"
            print("# falling back to entrez direct")
            print (cmd)
        else:
            for part in get_links(text):
                url2 = f"{url1}{part}"
                print (f"# {url2}")

                for link in get_links(read_url(url2)):
                    url3 = f"{url2}{link}"
                    print(f"wget -nc {url3}")








@plac.opt('db', "database")
@plac.opt('limit', "how many items to download")
@plac.opt('format_', "output format")
@plac.flg('quiet', "quiet mode")
def run(db='nuccore', format_="", quiet=False, limit=None, *data):
    """
    Writes code that fetches data from NCBI.
    """

    ids = []
    # Collect the ids
    for entry in data:
        if os.path.isfile(entry):
            elems = utils.read_lines(entry)
            ids.extend(elems)
        else:
            ids.append(entry)

    # A shortcut notation
    format_ = format_ if format_ else 'gbwithparts'

    # Writes the command
    write_command(ids, db=db, format_=format_)
