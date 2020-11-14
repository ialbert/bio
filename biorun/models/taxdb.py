import os, csv, json, gzip, sys, tarfile, re, codecs
import shutil
from itertools import islice
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils
from urllib import request
from biorun.libs import placlib as plac
from itertools import count

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

# Indentation character
INDENT = '  '

# Fields separator

SEP = ', '

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
    comm_dict = {}
    print(f"*** parsing: {name}")
    for index, elems in enumerate(stream):
        taxid, name, label = elems[0], elems[2], elems[6]
        if label == 'scientific name':
            # Store name, rank, common name, parent
            name_dict[taxid] = [name, "", "", ""]
        elif label == 'genbank common name':
            comm_dict[taxid] = name

    # Fill in common genbank names when exist
    for key, cname in comm_dict.items():
        if key in name_dict:
            sciname, rank, _1, _2 = name_dict[key]
            if sciname != cname:
                name_dict[key][2] = cname

    return name_dict


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
            name_dict[child][3] = parent

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
    name_dict = parse_names(fname, limit=limit)

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
        print("")
        table.commit()
        table.close()

    # Save the names into the database
    save_table(NAMES, name_dict)

    # Save the nodes.
    save_table(GRAPH, node_dict)

    print("*** saving the JSON model")
    json_path = os.path.join(utils.DATADIR, JSON_DB)

    # JSON will only have the graph and names.
    store = dict(NAMES=name_dict, GRAPH=node_dict, SYNONYMS={}, BACK={})
    fp = open(json_path, 'wt')
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
        print_node(node=node, depth=depth, names=names)
        visited.add(node)
        for neighbour in graph.get(node, []):
            dfs(visited, graph, neighbour, names, depth=depth + 1)


def print_node(node, names, depth):
    sep = SEP
    indent = INDENT * depth
    sciname, rank, cname, parent = names.get(node, ("MISSING", "NO RANK", "", ""))
    if cname and cname != sciname:
        print(f"{indent}{rank}{sep}{sciname} ({cname}){sep}{node}")
    else:
        print(f"{indent}{rank}{sep}{sciname}{sep}{node}")


def bfs(visited, graph, node, names):
    visited.add(node)
    queue.append(node)

    while queue:
        curr = queue.pop(0)

        print_node(curr, names=names)

        for nbr in graph.get(curr, []):
            if nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)


def backprop(node, names, collect=[]):
    if node in names:
        sciname, rank, cname, parent = names[node]
        if parent and parent != node:
            collect.append(parent)
            backprop(parent, names, collect)


def print_lineage(taxid, preload=False):
    names, graph = get_data(preload=preload)
    step = count(0)
    if taxid in names:
        collect = [taxid]
        backprop(taxid, names, collect=collect)
        for elem in reversed(collect):
            print_node(elem, names=names, depth=next(step))


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
    print(f"TaxDB nodes={node_size:,d} edges={edge_size:,d}")


def search_taxa(word, preload=False):
    names, graph = get_data(preload=preload)

    print(f"# searching taxonomy for: {word}")
    for taxid, name in search_names(word):
        print_node(taxid, names=names, depth=0)

def check_num(value):
    try:
        int(value)
        return True
    except ValueError as exc:
        return False

def query(taxid, preload=False):
    names, graph = get_data(preload=preload)

    isnum = check_num(taxid)
    if isnum and (taxid not in names):
        print(f"# taxid not found in database: {taxid}")
        sys.exit()

    if taxid in names:
        visited = set()
        dfs(visited, graph, taxid, names=names)
    else:
        search_taxa(taxid, preload=preload)

@plac.flg('build', "build a database from a taxdump")
@plac.flg('download', "download newest taxdump from NCBI")
@plac.flg('preload', "loads entire database in memory")
@plac.flg('lineage', "show the lineage for a taxon term", abbrev="L")
@plac.opt('indent', "the indentation string")
@plac.opt('sep', "separator string", abbrev="S")
@plac.opt('limit', "limit the number of entries", type=int)
@plac.flg('verbose', "verbose mode, prints more messages")
def run(limit=0, indent='   ', sep=', ', lineage=False, build=False, download=False, preload=False, verbose=False, *words):
    global SEP, INDENT

    limit = limit or None

    # Recognize string encodings: \t etc.
    INDENT = codecs.decode(indent, 'unicode_escape')
    SEP = codecs.decode(sep, 'unicode_escape')

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if download:
        download_taxdump()
    elif build:
        build_database(limit=limit)
    elif lineage:
        for word in words:
            print_lineage(word, preload=preload)
    else:
        # No terms listed. Print database stats.
        if not words:
            print_stats(preload=preload)

        # Print the query
        for word in words:
            query(word, preload=preload)


if __name__ == '__main__':
    # Bony fish: 117565
    # Betacoronavirus: 694002
    # SARS-COV2: 2697049

    # Jawless vertebrates: 1476529
    # Vertebrata: 7742

    plac.call(run)
