"""
Represents genome assemblies at NCBI
"""
from biorun.libs import placlib as plac
from biorun import const, utils, objects
import os, csv, json
from urllib.request import urljoin

join = os.path.join

ASSEMBLY_URL = "https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
ASSEMBLY_FILE_NAME = "assembly_summary_genbank.txt"

# Create the full paths
ASSEMBLY_FILE_NAME = join(utils.DATADIR, ASSEMBLY_FILE_NAME)


def resolve_fname(name, suffix='genomic.gbff.gz'):

    fname = os.path.join(utils.DATADIR, f"{name}_{suffix}")

    return fname


def fetch_data(url, acc, update=False):

    suffix = 'genomic.gbff.gz'

    end = f"{os.path.basename(url)}_{suffix}"
    url = f"{url}/{end}"
    dest = resolve_fname(acc, suffix=suffix)

    # File already exists
    if os.path.exists(dest) and not update:
        print(f"*** {acc} found at {dest}")
    else:
        utils.download(url=url, dest_name=dest)

    return


def search_file(fname=ASSEMBLY_FILE_NAME, accs=None, update=False):
    """
    Parse and search file
    """

    accs = [] if accs is None else accs

    if not os.path.isfile(fname):
        utils.error("assembly information not found (download first)")

    print(f"*** parsing {fname}")
    stream = open(fname)
    stream = filter(lambda x: x[0] != '#', stream)
    stream = csv.reader(stream, delimiter='\t')

    visited = set()
    accs = set(accs)
    nhits = 0

    def hit(a, r, g):
        hits = []
        # Iterate and collect hits
        for ac in accs:
            val = ac.strip()
            if (val == a or val == r or val == g) and ac not in visited:
                hits.append(val)
        return hits

    for line in stream:

        # Get a row, each keyed by respective column name.
        row = dict(zip(const.GENOMES_HEADER, line))

        acc = row['assembly_accession']
        root = acc.split('.')[0]
        gcf = row['gbrs_paired_asm']
        path = row['ftp_path']

        # Check to see if this assembly is in the query list
        matches = hit(acc, root, gcf)

        # Fetch the data if found
        if matches:
            fetch_data(url=path, acc=acc, update=update)
            visited.update([acc, root, gcf])
            nhits += 1

        # Stop reading through file once all queries have been found.
        if nhits == len(accs):
            break

    # Notify if some are not found
    missing = accs.difference(visited)

    if missing:
        print(f'*** {",".join(missing)} not found, skipped ...')


def update_assembly(url=ASSEMBLY_URL, dest_name=ASSEMBLY_FILE_NAME):
    """
    Downloads the assembly data.
    """
    utils.download(url=url, dest_name=dest_name)


@plac.pos("accs", "data names")
@plac.flg('update', "Update the genome assembly file")
@plac.flg('genome', "Execute", abbrev='x')
@plac.flg('download', "downloads the database from the remote site", abbrev='d')
@plac.flg('verbose', "verbose mode, prints more messages")
def run(update=False, download=False, genome=False, verbose=False, *accs):
    """
    Downloads genome assembly data from NCBI.
    """
    if download:
        update_assembly()
        return

    # Search for accession numbers and download assembly if any found.
    search_file(accs=accs, update=update)


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
