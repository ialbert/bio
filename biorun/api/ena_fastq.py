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
@plac.opt("out", help="optional output to the file", abbrev='o')
@plac.opt("limit", help="how many reads to download", abbrev='l', type=int)
def run(srr='ERR12058121', limit=10, out=''):

    # Sets the logger lever.
    logger.setLevel("INFO")

    # Set the prefix
    out = out or srr

    # The name of the directory
    dname = os.path.dirname(out)

    # Make directory dname
    if dname and not os.path.isdir(dname):
        logger.info(f"Creating directory: {dname}")
        os.makedirs(dname, exist_ok=True)

    # Is the limit set
    all = not limit

    # Obtain the metadata.
    meta = get_metadata(srr=srr)

    # Bad metadata may make this fail
    try:
        if all:
            n = int(meta[0].get('read_count', 0))
            b = meta[0].get('fastq_bytes', "0;0").split(";")
            b = map(float, b)
            b = map(lambda x: float(x)/(1024**3), b)
            b = ", ".join(map(lambda x: f"{x:.1f} GB", b))
            logger.info(f"Downloading {b} GB with {n:,} reads for {srr}")
        else:
            logger.info(f"Downloading {limit:,} reads for {srr}")
    except Exception as exc:
        logger.warning("metadata parsing problem, might still work ...")

    # Create URLs from metadata.
    urls = meta[0]['fastq_ftp'].split(';')
    urls = map(lambda x: f"https://{x}", urls)
    urls = list(urls)

    # Iterate over urls and downlad each file.
    for idx, url in enumerate(urls):

        fpath = f"{out}_{idx + 1}.fastq.gz"

        if all:
            # Download all reads
            pname = CMD.split()[0]
            exit_code = os.system(f"command -v {pname} > /dev/null 2>&1")
            if exit_code != 0:
                utils.error(f"Unable to run: {pname}", stop=False)
                utils.error(f"Installation: micromamba install aria2c")

            # Form the download command.
            cmd = f"{CMD} -o {fpath} {url}"
            logger.info(f"Running: {cmd}")
            #sys.stderr.flush()

            #continue

            # Run the download command.
            exit_code = os.system(cmd)
            if exit_code != 0:
                utils.error(f"Error when running: {cmd}")

        else:
            # Stream to a file


            # logger.info(f"Downloading {url}")

            # Open stream to remote gzipped files.
            stream = utils.get_gz_lines(url, limit=limit*4)

            # Open local gzip file.
            fp = gzip.open(fpath, mode='wb')

            logger.info(f"Saving to {fpath}")

            for line in stream:
                line = line.encode("utf-8")
                fp.write(line)


if __name__ == '__main__':
    plac.call(run)


