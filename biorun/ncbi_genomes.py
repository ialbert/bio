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
SQLITE_DB = "genome.db"
JSON_DB = "genome.json"

# Create the full paths
ASSEMBLY_FILE_NAME = join(utils.DATADIR, ASSEMBLY_FILE_NAME)

SQLITE_DB = join(utils.DATADIR, SQLITE_DB)
JSON_DB = join(utils.DATADIR, JSON_DB)

# Table names
ACC = "ACCESSIONS"

GCF = "GCF"

# Relate accession to asm_name
ASM = 'ASM'


def resolve_fname(name, suffix='genomic.gbff.gz'):

    fname = os.path.join(utils.DATADIR, f"{name}_{suffix}")

    return fname


def get_data(preload=False):

    if preload:
        if not os.path.isfile(JSON_DB):
            utils.error(f"File not found (you must build it first): {JSON_DB}")
        store = json.load(open(JSON_DB))
        accs = store[ACC]
        gca2gcf = store[GCF]
        asm = store[ASM]

    else:
        accs = utils.open_db(ACC, fname=SQLITE_DB)
        gca2gcf = utils.open_db(GCF, fname=SQLITE_DB)
        asm = utils.open_db(ASM, fname=SQLITE_DB)

    return accs, gca2gcf, asm


def parse_file(fname=ASSEMBLY_FILE_NAME):
    """
    Parse the file
    """
    if not os.path.isfile(fname):
        utils.error("assembly information not found (download first)")

    print(f"*** parsing {fname}")
    stream = open(fname)
    stream = filter(lambda x: x[0] != '#', stream)
    stream = csv.reader(stream, delimiter='\t')

    accs = dict()
    # Dict keyed by tax ids and point to a list of accession numbers.
    gca2gcf = dict()
    asm_names = dict()

    for line in stream:

        # Get a row, each keyed by respective column name.
        row = dict(zip(const.GENOMES_HEADER, line))

        acc = row['assembly_accession']
        tid = row['taxid']
        gcf = row['gbrs_paired_asm']
        asm = row['asm_name']
        path = row['ftp_path']

        accs[acc] = dict(tid=tid, path=path, gcf=gcf, asm=asm)

        gca2gcf[gcf] = acc
        asm_names[asm] = acc

    return accs, gca2gcf, asm_names


def build_database(fname=ASSEMBLY_FILE_NAME):

    # Check the source file exists
    if not os.path.isfile(fname):
        utils.error("assembly information not found (download first)")

    # Parse accession numbers and
    accs_dict, gca2gcf, asm_names = parse_file(fname=fname)

    def save(table, vals):
        utils.save_table(table, vals, fname=SQLITE_DB)

    save(ACC, accs_dict)

    save(GCF, gca2gcf)

    save(ASM, asm_names)

    print("*** saving the JSON model")

    # JSON will only have the graph and names.
    store = dict(ACC=accs_dict, GCF=gca2gcf, ASM=asm_names)
    fp = open(JSON_DB, 'wt')
    json.dump(store, fp, indent=4)
    fp.close()

    return accs_dict, gca2gcf, asm_names


def update_assembly(url=ASSEMBLY_URL, dest_name=ASSEMBLY_FILE_NAME):
    """
    Downloads the assembly data.
    """
    utils.download(url=url, dest_name=dest_name)


def fetch_data(accs, adict, gcfdict, asm):
    """
    Fetch accession data
    """

    # Build first if the database does not exist.
    suffix = 'genomic.gbff.gz'

    for acc in accs:
        # Search for this number in the GCF_, asm and accession dict
        gca = gcfdict.get(acc)
        asm_name = asm.get(acc)
        data = adict.get(acc) or adict.get(gca) or adict.get(asm_name)

        if not data:
            print(f"*** {acc} not found, skipping ...")
            continue

        url = data['path']
        end = f"{os.path.basename(url)}_{suffix}"
        url = f"{url}/{end}"

        # Where tp store the file
        dest = resolve_fname(acc, suffix=suffix)

        #if os.path.exists(dest):
        # Get the destination files
        # Download file into
        utils.download(url=url, dest_name=dest)

    return


@plac.pos("accs", "data names")
@plac.flg('update', "obtain the latest taxdump from NCBI")
@plac.flg('build', "Build the database")
@plac.flg('genome', "Execute", abbrev='x')
@plac.flg('download', "downloads the database from the remote site", abbrev='d')
@plac.flg('verbose', "verbose mode, prints more messages")
def run(update=False, download=False, build=False, genome=False, verbose=False, *accs):
    """
    Downloads genome assembly data from NCBI.
    """
    if download:
        update_assembly()
        return

    if build:
        # If the database does not exist, then
        build_database()
        return

    adict, gcf2gca, asm = get_data()

    fetch_data(accs=accs, adict=adict, gcfdict=gcf2gca, asm=asm)


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
