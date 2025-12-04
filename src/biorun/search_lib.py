"""
Bio search : search for information using accession numbers
"""
import re
import requests
import sys
import json
import csv
from biorun import logger, utils_lib

# Initialize the logger.
logger = logger.get_logger("main")

# SRR numbers: SRR5260547
SRR_PATT = re.compile(r'(ERR|SRR|DRR|SRP|ERP)\d+')
SRR_TYPE = 'srr'

# Bioproject numbers: PRJNA374918
PRJ_PATT = re.compile(r'PRJ([A-Z])+\d+')
PRJ_TYPE = 'bioproject'

# Genbank accessions: NC_045512
GBK_PATT = re.compile(r'(?P<letters>[a-zA-Z]+)(?P<under>_?)(?P<digits>\d+)(\.(?P<version>\d+))?')
GBK_TYPE = 'genbank'

# Ensemble ids: ENSG00000157764
ENS_PATT = re.compile(r'ENS[A-Z]+\d+')
ENS_TYPE = 'ensembl'

# A pattern that matches words followed by a colon and the literal wordgene
GENE_PATT = re.compile(r'(?P<word>\w+):gene', flags=re.IGNORECASE)
GENE_TYPE = 'gene'

# GEO accessions: GSM123456, GSE123456, GPL123456, GDS123456
GEO_PATT = re.compile(r'(GSM|GSE|GPL|GDS)\d+')
GEO_TYPE = 'geo'

# NCBI assembly ids: GCF_000001405.39, GCA_000001405.15
NCBI_ASM_PATT = re.compile(r'(GCA|GCF)_\d+')
NCBI_ASM_TYPE = 'assembly'

GENERIC_TYPE = 'generic'

# Patterns and types.
PATTERNS = [
    (NCBI_ASM_TYPE, NCBI_ASM_PATT),
    (GEO_TYPE, GEO_PATT),
    (GENE_TYPE, GENE_PATT),
    (SRR_TYPE, SRR_PATT),
    (PRJ_TYPE, PRJ_PATT),
    (ENS_TYPE, ENS_PATT),
    (GBK_TYPE, GBK_PATT),
]

# ENA PortalAPI urls.
# Documentation at https://www.ebi.ac.uk/ena/portal/api
ENA_PORTAL_API = "https://www.ebi.ac.uk/ena/portal/api"
ENA_PORTAL_FIELDS = f"{ENA_PORTAL_API}/returnFields"
ENA_PORTAL_REPORT = f"{ENA_PORTAL_API}/filereport"

# Default fields for ENA report.
ENA_DEFAULT_FIELDS = ['run_accession', 'sample_accession', 'sample_alias', 'sample_description', 
    'country', 'scientific_name',  'fastq_bytes',
    'library_name', 'library_strategy', 'library_source', 'library_layout', 'instrument_platform', 
    'instrument_model', 'study_title', 'first_public', 'base_count', 'read_count', 'fastq_ftp', ]

def get_valid_fields( url=ENA_PORTAL_FIELDS, db='ena'):
    """
    Returns all valid ENA fields for a database
    """
    params = dict(dataPortal=db, result='read_run')

    r = requests.get(url, params=params)

    if not r.ok:
        logger.error(f"Connection error: {r.status_code}")
        #r.raise_for_status()
        sys.exit(1)

    stream = r.iter_lines(decode_unicode=True)
    reader = csv.DictReader(stream, delimiter="\t") 
    return reader

def parse_fields(url=ENA_PORTAL_FIELDS):
    reader = get_valid_fields(url=url)
    fields = [r.get('columnId') for r in reader]
    fields = filter(lambda x: x, fields)
    fields = list(fields)
    return fields

def ena_srr_report(acc, all=False):
    """
    Get the ENA SRR report for an accession.
    """
  
    fields = parse_fields() if all else ENA_DEFAULT_FIELDS
    fields = ",".join(fields)
    params = dict(accession=acc, fields=fields, result='read_run')
    r = requests.get(ENA_PORTAL_REPORT, params=params)
    stream = r.iter_lines(decode_unicode=True)
    reader = csv.DictReader(stream, delimiter="\t")

    def add_field(entry):
        try:
            fastq_ftp = entry.get('fastq_ftp', ';')
            entry['fastq_ftp'] = fastq_ftp.split(";")

            fastq_bytes = entry.get('fastq_bytes', ';').split(";")
            fastq_bytes = map(lambda x: round(float(x)/1024/1024/1024,1), fastq_bytes)
            fastq_bytes = map(lambda x: f"{x} GB", fastq_bytes)
            fastq_bytes = list(fastq_bytes)
            fastq_bytes = ", ".join(fastq_bytes)
            entry['fastq_size'] = fastq_bytes
        except Exception as exc:
            logger.error(f"Error adding fields: {exc}")
        return entry

    reader = map(add_field, reader)
    return reader

# ENA ID lookup URL
ENA_ID_LOOKUP_URL = f"https://rest.ensembl.org/lookup/id"
def ena_lookup(acc, url=ENA_ID_LOOKUP_URL):
    """
    Connect to the lookup service.
    """
    conn = f"{url}/{acc}"
    
    try:    
        r = requests.get(conn, headers={"Content-Type": "application/json"})
    except Exception as exc:
        logger.error(f"Error: {exc}")
        return None

    if not r.ok:
        logger.error(f"Connection error: {r.status_code}")
        #r.raise_for_status()
        sys.exit(1)

    data = r.json()
    
    return data

def get_type(acc):
    """
    Get the type of the accession.
    """
    for dtype, pattern in PATTERNS:
        if pattern.match(acc):
            return dtype
    return GENERIC_TYPE

def as_json(data, indent=4):
    try:
        text = json.dumps(data, indent=indent)
    except Exception as exc:
        print(f"# JSON decoding error: {exc}")
        data = dict(error=f"{exc}")
        text = json.dumps(data, indent=indent)
       
    return text

PATT_FUNC = {
    NCBI_ASM_TYPE: None,
    GEO_TYPE: None,
    GENE_TYPE: None,
    SRR_TYPE: None,
    PRJ_TYPE: None,
    ENS_TYPE: ena_lookup,
    GBK_TYPE: None,
}

def run_single(acc):

    dtype = get_type(acc)
    func = PATT_FUNC.get(dtype, None)

    if not func:
        logger.error(f"# Can't recognize the accession: {acc}")
        return {}

    data = {}
    #data = func(acc)

    #print (data)

    return data


def run_all(accs):
    results = map(run_single, accs)
    results = list(results)
    return results

if __name__ == '__main__':
    accs = [
            "NP_001191", 
            "SRR5260547", 
            "PRJNA374918", 
            "HAD3:gene",
            "ENSG00000157764", 
            "GSM123456", 
            "GSE123456",
            "GCF_000001405.39",
            "GCA_000001405",
            "NP_001191.1", 
            "ecoli"
            ]
    
    #result = run_all(accs)

    #fields = get_valid_fields()
    #print(fields)
    reader = ena_srr_report("SRR5260547")

    #utils_lib.as_csv(reader)
    utils_lib.as_json(reader)
    


