# Needs:
# 
# pip install sqlitedict
#

import os, csv, json, gzip, sys, tarfile
import shutil
from itertools import islice
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils
from urllib import request
import plac

HOME = os.path.expanduser("~/.taxonkit")

join = os.path.join

nodes_fname = join(HOME, "nodes.dmp")
names_fname = join(HOME, "names.dmp")

# Create the thing here.

GRAPH, BACKLINK, NAMES = "GRAPH", "BACKLINK", "NAMES"

LIMIT = None

CHUNK = 10000

# Taxonomy database
TAXDB_URL = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

logger = utils.logger

def download_taxdump(url=TAXDB_URL, fname="taxdump.tar.gz"):
    """
    Downloads taxdump file.
    """
    print(f"*** downloading taxdump: {url}")
    path = os.path.join(utils.DATADIR, fname)

    # Download the data.
    src = request.urlopen(url)
    dest = open(path, 'wb')
    shutil.copyfileobj(src, dest)

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

    print("*** parsing names")
    for index, elems in enumerate(name_stream):
        taxid, name, label = elems[0], elems[2], elems[6]
        if label == 'scientific name':
            name_dict[taxid] = [name, ""]

    print("*** parsing nodes")
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
    graph.close()

    print ("")

def open_db(table, fname="taxonomy.db", flag='c'):
    """
    Opens a connection to a data table.
    """
    path = os.path.join(utils.DATADIR, fname)
    conn = SqliteDict(path, tablename=table, flag=flag, encode=json.dumps, decode=json.loads)
    return conn


queue = list()


def dfs(visited, graph, node, names, depth=0):
    if node not in visited:
        sciname, rank = names.get(node, ("MISSING", "NO RANK"))
        indent = "   " * depth
        print(f"{indent}{rank}, {sciname} ({node})")
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


def query(taxid):

    names = open_db(NAMES)
    graph = open_db(GRAPH)

    if taxid in names:
        visited = set()
        dfs(visited, graph, taxid, names=names)

    else:
        print(f"Invalid taxid: {taxid}")


@plac.flg('build', "download and build a new database")
@plac.flg('verbose', "verbose mode, prints more messages")
def run(build=False, verbose=False, *words):

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if build:
        download_taxdump()
    else:
        for word in words:
            query(word)

if __name__ == '__main__':
    # Bony fish: 117565
    # Betacoronavirus: 694002
    # SARS-COV2: 2697049

    plac.call(run)
