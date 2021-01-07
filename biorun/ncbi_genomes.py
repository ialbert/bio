"""
Represents genome assemblies at NCBI
"""
from biorun.libs import placlib as plac
from biorun import const, utils, objects
import os, csv

join = os.path.join

ASSEMBLY_URL = "https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
ASSEMBLY_FILE_NAME = "assembly_summary_genbank.txt"


# Create the full paths
ASSEMBLY_FILE_NAME = join(utils.DATADIR, ASSEMBLY_FILE_NAME)


def parse_file(fname=ASSEMBLY_FILE_NAME):
    """
    Parse the file
    """
    if not os.path.isfile(fname):
        utils.error("assembly information not found (download first)")

    stream = open(fname)
    stream = filter(lambda x: x[0] != '#', stream)
    stream = csv.reader(stream, delimiter='\t')
    for line in stream:
        print (line)
        break

def update_assembly(url=ASSEMBLY_URL, dest_name=ASSEMBLY_FILE_NAME):
    """
    Downloads the assembly data.
    """
    utils.download(url=url, dest_name=dest_name)


@plac.flg('update', "obtain the latest taxdump from NCBI")
@plac.flg('download', "downloads the database from the remote site", abbrev='d')
@plac.flg('verbose', "verbose mode, prints more messages")
def run(update=False, download=False, verbose=False, *words):
    """
    Downloads genome assembly data from NCBI.
    """
    if download:
        update_assembly()

    parse_file()

#
# Assembly summary information here
#
# https://ftp.ncbi.nih.gov/genomes/README_assembly_summary.txt
#
if __name__ == '__main__':
    """
    
    Genome data handling:

    Falciparium
    
    https://www.ncbi.nlm.nih.gov/assembly/GCF_000002765.5
    
    bio --fetch GCF_000002765
   
    bio --fetch GCF_000002765.5
   
    bio --fetch GCF_000002765.5 --gff --protein 

    """
    plac.call(run)
