'''
A simplified reimplementation of Entrez Direct interfaces

http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4

'''
import requests
from biorun import utils
from biorun.libs import xmltodict
from urllib.parse import urlsplit, urlunsplit
import json, os, csv
from biorun import const

# Entrez URL settings.
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI Assembly
ASSEMBLY_URL = "https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
ASSEMBLY_FILE_NAME = "assembly_summary_genbank.txt"

# Create the full paths
ASSEMBLY_FILE_NAME = os.path.join(utils.DATADIR, ASSEMBLY_FILE_NAME)

# Logging function
logger = utils.logger

# The keys that are valid environment keys.
ENV_KEYS = ["webenv", "use_history", "query_key"]


def fetch_genbank(acc, dest_name):
    """
    Returns a genbank file.
    """
    try:
        db = 'nuccore'

        rettype, retmode = "gbwithparts", "text"

        params = dict(db=db, rettype=rettype, id=acc, retmode=retmode)

        utils.download(EFETCH_URL, params=params, dest_name=dest_name)

    except Exception as exc:
        utils.error(exc)


def efetch(db, env={}, **kwds):
    """
    Perform an efetch query.
    Returns a dictionary based data structure.
    """
    try:

        params = dict(db=db)

        # Fill in with environment.
        for key, value in env.items():
            if key in ENV_KEYS:
                params[key] = value

        params.update(kwds)
        r = requests.get(EFETCH_URL, params=params)

        logger.info(r.url)

        data = xmltodict.parse(r.text)

        # Roundtrip to get rid of OrderedDicts
        data = json.loads(json.dumps(data))

        return data

    except Exception as exc:
        logger.error(exc)


def esearch(db, **kwds):
    """
    Performs an esearch query.
    Returns a dictionary based data structure.
    """
    try:
        params = dict(db=db)
        params.update(kwds)

        r = requests.get(ESEARCH_URL, params=params)

        logger.info(r.url)

        data = xmltodict.parse(r.text)

        # Roundtrip to get rid of OrderedDicts
        data = json.loads(json.dumps(data))

        # Parse the web environment
        field = data['eSearchResult']

        env = dict(
            webenv=field['WebEnv'],
            count=field['Count'],
            query_key=field['QueryKey'],
            use_history="y"
        )

        return env

    except Exception as exc:
        logger.error(exc)


def download_file(url, dest):
    # Split the urls
    (scheme, loc, path, query, frag) = urlsplit(url)

    # Last member of the path is the prefix for each file.
    prefix = path.split("/")[-1]

    # Remote file name.
    file_name = f"{prefix}_genomic.gbff.gz"

    # Remote file path.
    file_path = f"{path}/{file_name}"

    # Genbank url.
    file_url = urlunsplit((scheme, loc, file_path, '', ''))

    # Download the file.
    utils.download(url=file_url, dest_name=dest)

    return


def fetch_genome(acc, dest, update=False, fname=ASSEMBLY_FILE_NAME, ):
    """
    Parse and search and assembly file for an accession number.
    """

    # Update assembly information if it is missing.
    if not os.path.isfile(fname):
        update = True

    # Force the update
    if update:
        logger.info("updating assembly information")
        download_assembly()

    # Read the file line by line.
    logger.info(f"*** parsing {fname}")
    stream = open(fname, 'rt', encoding='utf-8')
    stream = filter(lambda x: x[0] != '#', stream)
    stream = csv.DictReader(stream, fieldnames=const.GENOME_ASSEMBLY_HEADER, delimiter='\t')

    # Read the through the file to to find the
    for row in stream:

        # Genbank version and root.
        gb_vers = row['assembly_accession']
        gb_base = gb_vers.split('.')[0]

        # Refseq version and root.
        rf_vers = row['gbrs_paired_asm']
        rf_base = rf_vers.split('.')[0]

        # The path to the file.
        url = row['ftp_path']

        # Found the match. Store by accession number.
        if acc in (gb_base, gb_vers, rf_base, rf_vers):
            download_file(url=url, dest=dest, acc=acc, update=update)
            return

    # If we go this far we have not found the data.
    print(f'*** accession not found: {acc}')


def download_assembly(url=ASSEMBLY_URL, dest_name=ASSEMBLY_FILE_NAME):
    """
    Downloads the assembly data.
    """
    utils.download(url=url, dest_name=dest_name)
