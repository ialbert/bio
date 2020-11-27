import os, csv, json, gzip, sys, tarfile, re, codecs
import shutil
from itertools import islice
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils
from urllib import request
from biorun.libs import placlib as plac
from itertools import count
from biorun.models import jsonrec
from biorun import storage

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

    if not os.path.isfile(fname):
        utils.error("taxdump file not found (download and build it first)")

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


def dfs(graph, node, names, depth=0, collect=[], visited=None):
    # Initialize the visited nodes once.
    visited = visited if visited else set()

    if node not in visited:
        text = node_formatter(node, names=names, depth=depth)
        print(text)
        collect.append((depth, node))
        visited.add(node)
        for nbr in graph.get(node, []):
            dfs(graph=graph, node=nbr, names=names, depth=depth + 1, collect=collect, visited=visited)


def get_values(node, names):
    # Avoiding code duplication everywhere
    sname, rank, cname, parent = names.get(node, ("MISSING", "NO RANK", "", ""))
    return sname, rank, cname, parent


def node_formatter(node, names, depth):
    """
    Creates a long form representation of a node.
    """
    sep = SEP
    indent = INDENT * depth
    sname, rank, cname, parent = get_values(node, names)

    # Decide what to do with common names.
    if cname and cname != sname:
        text = f"{indent}{rank}{sep}{sname} ({cname}){sep}{node}"
    else:
        text = f"{indent}{rank}{sep}{sname}{sep}{node}"

    return text


def backprop(node, names, collect=[]):
    if node in names:
        sciname, rank, cname, parent = names[node]
        if parent and parent != node:
            collect.append(parent)
            backprop(parent, names, collect)


def print_lineage(taxid, names, flat=1):
    step = count(0)
    if taxid in names:
        collect = [taxid]
        backprop(taxid, names, collect=collect)

        collect = collect[:-1]
        collect = reversed(collect)

        if flat:
            output = []
            for node in collect:
                sname, rank, cname, parent = get_values(node, names)
                output.append(sname)

            result = [taxid, ";".join(output)]
            print("\t".join(result))

        else:
            for node in collect:
                text = node_formatter(node, names=names, depth=next(step))
                print(text)


def get_data(preload=False):
    if preload:
        if not os.path.isfile(JSON_DB):
            utils.error(f"taxonomy file not found (you must build it first): {JSON_DB}")
        store = json.load(open(JSON_DB))
        names = store[NAMES]
        graph = store[GRAPH]
    else:
        names = open_db(NAMES)
        graph = open_db(GRAPH)

    return names, graph


def print_stats(names, graph):
    node_size, graph_size = len(names), len(graph)
    print(f"TaxDB: nodes={node_size:,d} parents={graph_size:,d}")


def search_taxa(word, preload=False):
    names, graph = get_data(preload=preload)

    word = codecs.decode(word, 'unicode_escape')

    print(f"# searching taxonomy for: {word}")
    for taxid, name in search_names(word):
        text = node_formatter(taxid, names=names, depth=0)
        print(text)


def check_num(value):
    try:
        int(value)
        return True
    except ValueError as exc:
        return False


def print_database(names, graph):
    for name in names:
        text = node_formatter(name, names=names, depth=0)
        print(text)


def query(taxid, names, graph):
    """
    Prints the descendants of node
    """
    isnum = check_num(taxid)
    if isnum and (taxid not in names):
        print(f"# taxid not found in database: {taxid}")
        sys.exit()

    if taxid in names:
        collect = []
        dfs(graph, taxid, names=names, collect=collect)


    else:
        search_taxa(taxid)


@plac.pos("words", "taxids or search queries")
@plac.flg('build', "build a database from a taxdump")
@plac.flg('download', "download newest taxdump from NCBI")
@plac.flg('preload', "loads entire database in memory")
@plac.flg('list_', "lists database content")
@plac.flg('flat', "flattened output")
@plac.flg('lineage', "show the lineage for a taxon term", abbrev="L")
@plac.opt('indent', "the indentation string")
@plac.opt('sep', "separator string", abbrev="S")
@plac.opt('limit', "limit the number of entries", type=int, abbrev='T')
@plac.flg('verbose', "verbose mode, prints more messages")

def run(limit=0, list_=False, flat=False, indent='   ', sep=', ', lineage=False, build=False, download=False,
        preload=False,
        verbose=False, *words):
    global SEP, INDENT

    limit = limit or None

    # Recognize string encodings: \t etc.
    INDENT = codecs.decode(indent, 'unicode_escape')
    SEP = codecs.decode(sep, 'unicode_escape')

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Access the database.
    names, graph = get_data(preload=preload)

    if list_:
        print_database(names=names, graph=graph)
        sys.exit()

    if download:
        download_taxdump()

    if build:
        build_database(limit=limit)

    terms = []
    # Attempts to fetch data if possible.
    for word in words:
        json = storage.get_json(word)
        doubles = [jsonrec.find_taxid(rec) for rec in json] if json else [[]]
        taxids = [elem for sublist in doubles for elem in sublist]
        if taxids:
            terms.extend(taxids)
        else:
            terms.append(word)

    for word in terms:

        if lineage:
            print_lineage(word, names=names, flat=flat)
        else:
            query(word, names=names, graph=graph)

    # No terms listed. Print database stats.
    if not terms:
        print_stats(names=names, graph=graph)


if __name__ == '__main__':
    # Bony fish: 117565
    # Betacoronavirus: 694002
    # SARS-COV2: 2697049

    # Jawless vertebrates: 1476529
    # Vertebrata: 7742

    plac.call(run)
