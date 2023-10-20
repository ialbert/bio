"""
Connects to the ENA (European Nucleotide Archive) to download FASTQ files.
"""

import json, sys, gzip, os
from biorun import utils
from biorun.libs import placlib as plac

# Use the package logger
logger = utils.logger

# ENA API points
ENA_API = "https://www.ebi.ac.uk/ena/portal/api"
ENA_FIELDS = f"{ENA_API}/returnFields"
ENA_REPORT = f"{ENA_API}/filereport"

# The download command with aria2c
CMD = 'aria2c -x 5 -c --summary-interval 10'

# Fetch metadata
def get_metadata(srr):

    # The SRR number
    # logger.info(f"Searching Ensembl for {srr}")

    # The URL to fetch
    url = ENA_REPORT

    # The fields to fetch
    fields = [ "run_accession", 'fastq_ftp', 'fastq_md5', 'fastq_bytes',  'read_count' ]

    # Form the fields string
    fields = ",".join(fields)

    # The parameters for the query
    params = dict(
        accession=srr, fields=fields,
        format='json', result='read_run',
    )

    # Download the metadata
    resp = utils.get_url(url, params=params)

    # Decode the JSON
    try:
        resp = json.loads(resp)
    except Exception as exc:
        utils.error(f"JSON decoding error: {exc}")

    return resp


def test_metadata(srr='ERR12058121'):
    """
    Testing the metadata fetch.
    """
    meta = get_metadata(srr=srr)

    print(json.dumps(meta, indent=4))


@plac.pos("srr", help="the srr numbers", )
def run(srr, **kwargs):

    # Obtain the metadata.
    meta = get_metadata(srr=srr)

    try:
        # Create URLs from metadata.
        urls = meta[0]['fastq_ftp'].split(';')
        urls = map(lambda x: f"https://{x}", urls)
        urls = list(urls)

        # Display the urls for the data.
        for url in urls:
            print(url)
    except Exception as exc:
        utils.error(f"Metadata parsing error: {exc}")


if __name__ == '__main__':
    plac.call(run)


