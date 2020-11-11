
import os, csv, json, gzip, sys, tarfile
import shutil
from itertools import islice
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils
from urllib import request
import plac

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

GRAPH, BACKLINK, NAMES = "GRAPH", "BACKLINK", "NAMES"

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



def build_database(url=TAXDB_URL, fname=TAXDB_NAME):
    """
    Downloads taxdump file.
    """
    print(f"*** building database from: {fname}")
    path = os.path.join(utils.DATADIR, fname)

    if not os.path.isfile(path):
        utils.error(f"no taxdump file found")

    # The taxdump file.
    tar = tarfile.open(path, "r:gz")

    def get_stream(name):
        mem = tar.getmember(name)
        stream = tar.extractfile(mem)
        stream = map(lambda x: x.decode("ascii"), stream)
        stream = islice(stream, LIMIT)
        return stream

    name_stream = get_stream("names.dmp")
    name_stream = csv.reader(name_stream, delimiter="\t")

    node_stream = get_stream("nodes.dmp")
    node_stream = csv.reader(node_stream, delimiter="\t")

    node_dict = {}
    name_dict = {}
    back_dict = {}

    print("*** parsing: names.dmp")
    for index, elems in enumerate(name_stream):
        taxid, name, label = elems[0], elems[2], elems[6]
        if label == 'scientific name':
            name_dict[taxid] = [name, ""]

    print("*** parsing: nodes.dmp")
    for elems in node_stream:
        child, parent, rank = elems[0], elems[2], elems[4]
        back_dict[child] = parent
        node_dict.setdefault(parent, []).append(child)
        if child in name_dict:
            name_dict[child][1] = rank

    # Save the names into the database
    nsize = len(name_dict)
    names = open_db(table=NAMES, flag='w')
    for index, (key, value) in enumerate(name_dict.items()):
        names[key] = value
        if index % CHUNK == 0:
            perc = round(index / nsize * 100)
            print(f"*** saving {nsize:,} names ({perc:.0f}%)", end="\r")
            names.commit()
    names.commit()
    names.close()

    print("")

    # Save the nodes.
    gsize = len(node_dict)
    graph = open_db(table=GRAPH, flag='w')
    for index, (key, value) in enumerate(node_dict.items()):
        graph[key] = value
        if index % CHUNK == 0:
            perc = round(index / gsize * 100)
            print(f"*** saving {gsize:,} nodes ({perc:.0f}%)", end="\r")
            graph.commit()
    graph.commit()
    graph.close()

    print ("")

    print ("*** saving the JSON model")
    json_path = os.path.join(utils.DATADIR, JSON_DB)
    store = dict(NAMES=name_dict, GRAPH=node_dict)
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
        print(f"{indent}{rank}, {sciname}, taxid={node}")
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


def query(taxid, mode=False):

    if mode:
        store = json.load(open(JSON_DB))
        names = store[NAMES]
        graph = store[GRAPH]
    else:
        names = open_db(NAMES)
        graph = open_db(GRAPH)


    if taxid in names:
        visited = set()
        dfs(visited, graph, taxid, names=names)

    else:
        print(f"Invalid taxid: {taxid}")


@plac.flg('build', "build a database from a taxdump")
@plac.flg('download', "download newest taxdump from NCBI")
@plac.flg('preload', "loads entire database in memory")
@plac.flg('verbose', "verbose mode, prints more messages")
def run(build=False, download=False, preload=False, verbose=False, *words):

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if download:
        download_taxdump()

    if build:
        build_database()

    for word in words:
        query(word, mode=preload)

if __name__ == '__main__':
    # Bony fish: 117565
    # Betacoronavirus: 694002
    # SARS-COV2: 2697049

    # Jawless vertebrates: 1476529
    # Vertebrata: 7742

    plac.call(run)
