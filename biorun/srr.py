'''
Module that deals with SRR numbers.
'''
from biorun.libs import placlib as plac
import re, json, sys, gzip, os, subprocess
from biorun import utils
from itertools import *
import requests
from urllib.parse import urlparse

# SRR numbers: SRR5260547
SRR = re.compile(r'(ERR|SRR|DRR|SRP|ERP)\d+')

logger = utils.logger

# ENA API points
ENA_API = "https://www.ebi.ac.uk/ena/portal/api"
ENA_FIELDS = f"{ENA_API}/returnFields"
ENA_REPORT = f"{ENA_API}/filereport"

# The download command with aria2c
CMD = 'aria2c -x 5 -c --summary-interval 10'

# Match SRR numbers
def match_srr(text):
    """
    Pattern for SRR numbers.
    """
    return bool(SRR.search(text))

def sumbit_request(url, params={}):

    try:
        r = requests.get(url, params=params)
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        return r.text

    except Exception as exc:
        utils.error(f"Error when searching for {srr}")


def read_lines(url, N=100):

    try:
        r = requests.get(url, stream=True)

        r.raise_for_status()

        stream = gzip.open(r.raw, mode='rt', encoding='utf-8')
        stream = islice(stream, N)

        for line in stream:
            yield (line)

    except Exception as exc:
        print(f"# Error: {exc}")

# Fetch metadata
def get_metadata(srr):
    logger.info(f"Searching Ensembl for {srr}")

    url = ENA_REPORT

    fields = [
        'run_accession',
        "sample_accession",
        "sample_alias",
        "sample_description",
        'first_public',
        'country',
        'scientific_name',
        'fastq_bytes',
        'base_count',
        'read_count',
        'library_name',
        "library_strategy",
        "library_source",
        'library_layout',
        'instrument_platform',
        'instrument_model',
        'study_title',
        'fastq_ftp',
    ]

    fields = ",".join(fields)

    params = dict(
        accession=srr,
        fields=fields,
        format='json',
        result='read_run',
    )

    resp = sumbit_request(url, params=params)

    try:
        resp = json.loads(resp)
    except Exception as exc:
        utils.error(f"JSON decoding error: {exc}")

    return resp


def test_metadata(srr='ERR12058121'):

    meta = get_metadata(srr=srr)

    fastq_ftps = meta[0]['fastq_ftp'].split(';')

    print (fastq_ftps)


@plac.pos("srr", help="the srr numbers", )
@plac.opt("dname", help="the directory name", abbrev='d')
@plac.opt("prefix", help="optional prefix to the filenames", abbrev='p')
@plac.opt("n", help="how many reads to download", abbrev='n')
@plac.flg("all", help="download all reads", abbrev='a')
def run(srr='ERR12058121', dname='reads', n=12, prefix='', all=False):

    # Set the prefix
    prefix = prefix or srr

    logger.setLevel("INFO")

    # Make directory dname
    os.makedirs(dname, exist_ok=True)

    meta = get_metadata(srr=srr)
    urls = meta[0]['fastq_ftp'].split(';')
    urls = map(lambda x: f"https://{x}", urls)
    urls = list(urls)

    for idx, url in enumerate(urls):

        fname = f"{prefix}_{idx+1}.fastq.gz"
        fpath = os.path.join(dname, fname)

        if all:
            # Download all reads
            logger.info(f"Downloading all reads with aria2")
            pname = CMD.split()[0]
            exit_code = os.system(f"command -v {pname} > /dev/null 2>&1")
            if exit_code != 0:
                logger.error(f"Unable to run: {pname}")
                logger.error(f"Installation: micromamba install aria2c")
                sys.exit(1)

            cmd = f"{CMD} -o {fpath} {url}"
            logger.info(f"Running: {cmd}")
            sys.stderr.flush()
            exit_code = os.system(cmd)
            if exit_code != 0:
                logger.error(f"Error when running: {cmd}")
                sys.exit(1)

        else:
            # Stream to a file
            logger.info(f"Limit to {n} reads")

            logger.info(f"Downloading {url}")

            fp = gzip.open(fpath, mode='wb')
            stream = read_lines(url, N=4)

            logger.info(f"Saving to {fpath}")

            for line in stream:
                line = line.encode("utf-8")
                fp.write(line)



if __name__ == '__main__':
    plac.call(run)


