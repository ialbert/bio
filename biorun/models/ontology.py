import json, shutil, os, tarfile, re
from urllib import request
from itertools import islice

from biorun.libs import placlib as plac
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils


JSON_DB = "ontology.json"
SQLITE_DB = "ontology.sqlite"
GENE_URL = "http://purl.obolibrary.org/obo/go.obo"
SEQ_URL = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so-simple.obo"
ONOTO_NAME = "ontology.obo"

# Delimiter in ontology files, separating one term from another.
DELIM = '[Term]'

# GO or SO id patterns.
ID_PATT = r"[G|S]O:\d+"

GENE_ID = 'GO'
SEQ_ID = 'SO'

join = os.path.join

# Create the full paths
ONOTO_NAME = join(utils.DATADIR, ONOTO_NAME)
SQLITE_DB = join(utils.DATADIR, SQLITE_DB)
JSON_DB = join(utils.DATADIR, JSON_DB)

# Table names

# Node descriptions
TERM = "TERMS"

# Graph following the is_a relation
GRAPH = "GRAPH"

# Names mapped to a GO id.
NAMES = "NAMES"

CHILDREN = "CHILDREN"

logger = utils.logger


CHUNK = 25000


def download_ontology(url, fname=ONOTO_NAME, mode="wb"):
    """
    Downloads ontology file.
    """
    print(f"*** downloading ontology: {url}")

    # Download the data.
    src = request.urlopen(url)
    dest = open(fname, mode=mode)
    shutil.copyfileobj(src, dest)


def download_terms():
    download_ontology(url=GENE_URL)

    # Append Sequence ontology to existing file.
    download_ontology(url=SEQ_URL, mode="ab")
    return


def get_data(preload=False):
    if preload:
        if not os.path.isfile(JSON_DB):
            utils.error(f"ontology file not found (you must build it first): {JSON_DB}")
        store = json.load(open(JSON_DB))
        terms = store[TERM]
        nodes = store[GRAPH]
        names = store[NAMES]
        back = store[CHILDREN]
    else:
        terms = open_db(TERM)
        nodes = open_db(GRAPH)
        names = open_db(NAMES)
        back = open_db(CHILDREN)

    return terms, nodes, names, back


def open_db(table, fname=SQLITE_DB, flag='c'):
    """
    Opens a connection to a data table.
    """
    conn = SqliteDict(fname, tablename=table, flag=flag, encode=json.dumps, decode=json.loads)
    return conn


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


def match_id(item):

    match = re.search(ID_PATT, item)

    oid = match.group(0) if match else None

    return oid


def stop_term(elem):
    stop = (DELIM in elem) or ("[Typedef]" in elem)
    return stop


def parse_term(fname):
    """
    Parses the ontology file into a terms dictionary.
    """

    # The ontology file, with both sequence and gene info.
    stream = open(fname, mode="r")

    terms = {}
    names = {}
    nodes = {}
    back_prop = {}
    uid, parent, name, definition = None, None, None, None

    print(f"*** parsing: {fname}")
    for elems in stream:

        if stop_term(elems):
            # Reset objs for next term.
            uid, parent, name, definition = None, None, None, None
            continue

        val = elems.split(":")

        if elems.startswith("id:"):
            uid = match_id(elems)

        if elems.startswith("is_a:"):
            parent = match_id(elems)

        if "name:" in elems:
            name = val[1].strip()

        if elems.startswith("def:"):
            definition = ":".join(val[1:]).strip()

        if uid:
            terms[uid] = [name, definition]
            names[name] = uid
            nodes.setdefault(parent, []).append(uid)
            back_prop[uid] = parent

    return terms, nodes, names, back_prop


def build_database():
    """
    Build the ontology database.
    """
    print(f"*** building database from: {ONOTO_NAME}")

    # Check the file.
    if not os.path.isfile(ONOTO_NAME):
        logger.info(f"No ontology file found, downloading...")
        download_terms()

    # Parse the terms from file
    terms, nodes, names, back_prop = parse_term(fname=ONOTO_NAME)

    # Save terms into the database
    save_table(TERM, terms)

    # Save graph into the database
    save_table(GRAPH, nodes)

    # Save names into the database
    save_table(NAMES, names)

    # Save child-parent into the database
    save_table(CHILDREN, back_prop)

    print("*** saving the JSON model")
    json_path = os.path.join(utils.DATADIR, JSON_DB)

    store = dict(TERMS=terms, GRAPH=nodes, NAMES=names, CHILDREN=back_prop)
    fp = open(json_path, 'wt')
    json.dump(store, fp, indent=4)
    fp.close()


def walk_tree(nodes, start, collect=None):

    collect = [] if collect is None else collect
    collect.append(start)
    parent = nodes.get(start)

    if parent:
        walk_tree(nodes=nodes, start=parent, collect=collect)


def print_tree(terms, tree=None):

    tree = [] if tree is None else tree
    tree = reversed(tree)

    # Print the definition and all children here.
    for idx, oid in enumerate(tree):
        vals = terms.get(oid)
        if vals:
            name, definition = vals
            pad = '\t' * idx
            print(f"{pad}{oid} {name}")


def printer(uids, terms):

    for uid in uids:
        name, definition = terms[uid]
        print(f"{uid}\t{name}")


def show_linage(start, terms, nodes):

    collect = []
    walk_tree(nodes=nodes, start=start, collect=collect)

    # Print the tree.
    print_tree(tree=collect, terms=terms)

    return


def search(query, terms, prefix=""):

    # Search for names containing query.

    for uid, vals in terms.items():
        name, definition = vals

        if prefix and not uid.startswith(prefix):
            continue

        # Print all terms containing this name.
        if (query in name) or (query in uid):
            name, definition = terms[uid]
            print(f"{uid} {name}")


def perform_query(query, terms, nodes, names, back_prop, prefix="", linage=False):
    """
    Search database based on a name or the ontological id.
    """

    # The query is an exact match, print info.
    if names.get(query) or terms.get(query):
        # Get the GO or SO id
        uid = names.get(query) or query
        if prefix and not uid.startswith(prefix):
            return

        # Fetch the name and definition
        name, definition = terms[uid]
        print(f"{uid}\t{name}\t{definition}")
        if linage:
            show_linage(start=uid, terms=terms, nodes=back_prop)
            return

        # Print children
        children = set(nodes.get(uid, []))
        printer(uids=children, terms=terms)

        return

    # Search for names containing query.
    search(query=query, terms=terms, prefix=prefix)


@plac.pos('query', "Search database by ontological name or GO/SO ids.")
@plac.flg('build', "build a database of all gene and sequence ontology terms. ")
@plac.flg('preload', "loads entire database in memory")
@plac.flg('verbose', "verbose mode, prints more messages")
@plac.flg('linage', "show the ontological linage")
@plac.flg('download', "download newest ontological data")
@plac.flg('so', "Filter query for sequence ontology terms.")
@plac.flg('go', "Filter query for gene ontology terms.")
def run(query="", build=False, download=False, preload=False, so=False, go=False,
        linage=False, verbose=False):

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if download:
        download_terms()

    if build:
        build_database()

    terms, nodes, names, back_prop = get_data(preload=preload)
    query = query.strip()

    prefix = SEQ_ID if so else ''
    prefix = GENE_ID if go else prefix

    if query:
        perform_query(query=query,
                      linage=linage,
                      terms=terms,
                      prefix=prefix,
                      nodes=nodes,
                      back_prop=back_prop,
                      names=names)

