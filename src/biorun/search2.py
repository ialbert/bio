"""
Bio search : search for information using accession numbers
"""
import re
from biorun.api import ena

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
ASM_PATT = re.compile(r'(GCA|GCF)_\d+')
ASM_TYPE = 'assembly'

# Patterns and types.
PATTERNS = [
    (ASM_TYPE, ASM_PATT),
    (GEO_TYPE, GEO_PATT),
    (GENE_TYPE, GENE_PATT),
    (SRR_TYPE, SRR_PATT),
    (PRJ_TYPE, PRJ_PATT),
    (ENS_TYPE, ENS_PATT),
    (GBK_TYPE, GBK_PATT),
]

# ENA API urls.
ENA_API = "https://www.ebi.ac.uk/ena/portal/api"
ENA_FIELDS = f"{ENA_API}/returnFields"
ENA_REPORT = f"{ENA_API}/filereport"

def get_match(acc):
    for dtype, pattern in PATTERNS:
        if pattern.match(acc):
            return dtype
    return None

def run(acc):

    dtype = get_match(acc)

    if not dtype:
        print(f"# No match for: {acc}")
    else:
        print (acc, dtype)

    if dtype == ENS_TYPE:
        #ena.lookup(acc, url=ena.LOOKUP_URL)
        pass

    if dtype == ASM_TYPE:
        ena.lookup(acc, url=ena.ASSEMBLY_URL)

    print("-"* 10)


if __name__ == '__main__':
    accs = [
            "NP_001191", "SRR5260547", "PRJNA374918", "HAD3:gene",
            "ENSG00000157764", "GSM123456", "GSE123456",
            "GCF_000001405.39",
            "GCA_000001405",
            "NP_001191.1", "ecoli"
            ]
    for acc in accs:
        run(acc)


