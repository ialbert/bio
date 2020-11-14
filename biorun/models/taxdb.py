
import os, csv, json, gzip, sys, tarfile, re
import shutil
from itertools import islice
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils
from urllib import request
from biorun.libs import placlib as plac

JSON_DB = "taxdb.json"
SQLITE_DB = "taxdb.sqlite"
TAXDB_URL = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAXDB_NAME = "taxdump.tar.gz"

join = os.path.join

# Create the full paths
TAXDB_NAME = join(utils.DATADIR, TAXDB_NAME)
SQLITE_DB = join(utils.DATADIR, SQLITE_DB)
JSON_DB = join(utils.DATADIR, JSON_DB)


# Create the thing here.

GRAPH, BACK, NAMES, SYNON = "GRAPH", "BACKLINKS", "NAMES", "SYNONYMS"

LIMIT = None

CHUNK = 25000

logger = utils.logger


def download_taxdump(url=TAXDB_URL, fname=TAXDB_NAME):
    """
    Downloads taxdump file.
    """
    print(f"*** downloading taxdump: {url}")

    # Download the data.
    src = request.urlopen(url)
    dest = open(fname, 'wb')
    shutil.copyfileobj(src, dest)

def get_stream(tar, name, limit=None):
    mem = tar.getmember(name)
    stream = tar.extractfile(mem)
    stream = map(lambda x: x.decode("ascii"), stream)
    stream = islice(stream, limit)
    return stream

def search_names(word, fname=TAXDB_NAME, name="names.dmp", limit=None):
    """
    Parses the names.dmp component of the taxdump.
    """

    # The taxdump file.
    tar = tarfile.open(fname, "r:gz")

    stream = get_stream(tar=tar, name=name, limit=limit)
    stream = csv.reader(stream, delimiter="\t")

    patt = re.compile(word, re.IGNORECASE)

    for index, elems in enumerate(stream):
        taxid, name, label = elems[0], elems[2], elems[6]
        if label == 'scientific name' or label == 'equivalent name' or label == 'genbank common name':
            if patt.search(name):
                yield taxid, name

def parse_names(fname, name="names.dmp", limit=None):
    """
    Parses the names.dmp component of the taxdump.
    """

    # The taxdump file.
    tar = tarfile.open(fname, "r:gz")

    stream = get_stream(tar=tar, name=name, limit=limit)
    stream = csv.reader(stream, delimiter="\t")

    name_dict = {}
    syn_dict  = {}

    print(f"*** parsing: {name}")
    for index, elems in enumerate(stream):
        taxid, name, label = elems[0], elems[2], elems[6]
        if label == 'scientific name':
            name_dict[taxid] = [name, ""]
            syn_dict[name] = taxid
        elif label == 'equivalent name':
            syn_dict[name] = taxid

    return name_dict, syn_dict


def parse_nodes(fname, name_dict, name="nodes.dmp", limit=None):
    """
    Parses the names.dmp component of the taxdump.
    """

    # The taxdump file.
    tar = tarfile.open(fname, "r:gz")

    stream = get_stream(tar=tar, name=name, limit=limit)
    stream = csv.reader(stream, delimiter="\t")

    node_dict = {}
    back_dict = {}

    print("*** parsing: nodes.dmp")
    for elems in stream:
        child, parent, rank = elems[0], elems[2], elems[4]
        back_dict[child] = parent
        node_dict.setdefault(parent, []).append(child)
        if child in name_dict:
            name_dict[child][1] = rank

    return node_dict, back_dict

def build_database(fname=TAXDB_NAME, limit=None):
    """
    Downloads taxdump file.
    """
    print(f"*** building database from: {fname}")
    path = os.path.join(utils.DATADIR, fname)

    # Check the file.
    if not os.path.isfile(path):
        utils.error(f"no taxdump file found, run the --download flag")

    # Parse the names
    name_dict, syn_dict = parse_names(fname, limit=limit)

    # Parse the nodes.
    node_dict, back_dict = parse_nodes(fname, name_dict=name_dict, limit=limit)

    def save_table(name, obj):
        size = len(obj)
        table = open_db(table=name, flag='w')
        for index, (key, value) in enumerate(obj.items()):
            table[key] = value
            if index % CHUNK == 0:
                perc = round(index / size * 100)
                print(f"*** saving {name} with {size:,} elements ({perc:.0f}%)", end="\r")
                table.commit()
        print(f"*** saved {name} with {size:,} elements (100%)", end="\r")
        print ("")
        table.commit()
        table.close()

    # Save the names into the database
    save_table(NAMES, name_dict)

    # Save the nodes.
    save_table(GRAPH, node_dict)

    print ("*** saving the JSON model")
    json_path = os.path.join(utils.DATADIR, JSON_DB)

    # JSON will only have the graph and names.
    store = dict(NAMES=name_dict, GRAPH=node_dict, SYNONYMS={}, BACK={})
    fp = open(json_path,'wt')
    json.dump(store, fp, indent=4)
    fp.close()

def open_db(table, fname=SQLITE_DB, flag='c'):
    """
    Opens a connection to a data table.
    """
    conn = SqliteDict(fname, tablename=table, flag=flag, encode=json.dumps, decode=json.loads)
    return conn


queue = list()


def dfs(visited, graph, node, names, depth=0):
    if node not in visited:
        sciname, rank = names.get(node, ("MISSING", "NO RANK"))
        indent = "   " * depth
        print(f"{indent}{rank}, {sciname}, {node}")
        visited.add(node)
        for neighbour in graph.get(node, []):
            dfs(visited, graph, neighbour, names, depth=depth + 1)


def bfs(visited, graph, node, names, ):
    visited.add(node)
    queue.append(node)

    while queue:
        tgt = queue.pop(0)
        sciname, rank = names.get(tgt, ("MISSING", "NO RANK"))
        print(f"{tgt}, {rank}, {sciname}")

        for nbr in graph.get(tgt, []):
            if nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)


def get_data(preload=False):

    if preload:
        store = json.load(open(JSON_DB))
        names = store[NAMES]
        graph = store[GRAPH]
    else:
        names = open_db(NAMES)
        graph = open_db(GRAPH)

    return names, graph

def print_stats(preload=False):
    names, graph = get_data(preload=preload)
    node_size, edge_size = len(names), len(graph)
    print (f"TaxDB nodes={node_size:,d} edges={edge_size:,d}")

def query(taxid, preload=False):


    names, graph = get_data(preload=preload)

    if taxid in names:
        visited = set()
        dfs(visited, graph, taxid, names=names)

    else:
        print(f"# searching taxonomy for: {taxid}")
        for taxid, name in search_names(taxid):
            sciname, rank = names.get(taxid, ('', ''))
            if sciname != name:
                print(f"{taxid}\t{rank}\t{sciname} ({name})")
            else:
                print (f"{taxid}\t{rank}\t{sciname}")

@plac.flg('build', "build a database from a taxdump")
@plac.flg('download', "download newest taxdump from NCBI")
@plac.flg('preload', "loads entire database in memory")
@plac.opt('limit', "limit the number of entries", type=int)
@plac.flg('verbose', "verbose mode, prints more messages")
def run(limit=0, build=False, download=False, preload=False, verbose=False, *words):

    limit = limit or None

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if download:
        download_taxdump()

    if build:
        build_database(limit=limit)

    for word in words:
        query(word, preload=preload)

    if not words:
        print_stats(preload=preload)

if __name__ == '__main__':
    # Bony fish: 117565
    # Betacoronavirus: 694002
    # SARS-COV2: 2697049

    # Jawless vertebrates: 1476529
    # Vertebrata: 7742

    plac.call(run)
