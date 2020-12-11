import json, shutil, os, tarfile, re
import requests
from itertools import islice
from textwrap import wrap

from biorun.libs import placlib as plac
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils


JSON_DB = "ontology.json"
SQLITE_DB = "ontology.sqlite"

GO_URL = "http://purl.obolibrary.org/obo/go.obo"
SO_URL = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so-simple.obo"
ONTOLOGY_FILE = "ontology.obo"

# Delimiter in ontology files, separating one term from another.
DELIM = '[Term]'

# GO or SO id patterns.
ID_PATT = r"[G|S]O:\d+"

# Edge type pattern
EDGE_TYPE_PATT = r"relationship:\s+(.+)(\s+[S|G])"

GO_ID = 'GO'
SO_ID = 'SO'

join = os.path.join

# Create the full paths
ONTOLOGY_FILE = join(utils.DATADIR, ONTOLOGY_FILE)
SQLITE_DB = join(utils.DATADIR, SQLITE_DB)
JSON_DB = join(utils.DATADIR, JSON_DB)


INDENT = '  '

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


def download_ontology(url, fname=ONTOLOGY_FILE, mode="wt"):
    """
    Downloads ontology file.
    """
    print(f"*** downloading: {url}")

    # Stream data into a file.
    r = requests.get(url, stream=True)
    fp = open(fname, mode)
    for chunk in r.iter_content(chunk_size=2**12):
        fp.write(chunk.decode("ascii", "backslashreplace"))

def download_terms():

    # Download Gene Ontology (GO).
    download_ontology(url=GO_URL)

    # Append Sequence Ontology (SO).
    download_ontology(url=SO_URL, mode="at")
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


def stop_term(text):
    stop = (DELIM in text) or ("[Typedef]" in text)
    return stop


def edge_type(item):
    """
    Extract edge type from the element
    """

    match = re.search(EDGE_TYPE_PATT, item)

    etype = match.group(1) if match else None

    return etype


def update_nodes(nodes, back_prop, edges):

    if not edges:
        return

    def update(par, child, et):
        nodes.setdefault(par, []).append((child, et))
        back_prop.setdefault(child, []).append((par, et))
        return

    for parent, current, etype in edges:
        update(parent, current, etype)

    return


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
    uid, parent, name, definition, edges = None, None, None, None, []

    print(f"*** parsing: {fname}")
    for elems in stream:

        if stop_term(elems):
            # Reset objs for next term.
            uid, parent, name, definition, edges = None, None, None, None, []
            continue

        val = elems.split(":")

        if elems.startswith("id:"):
            uid = match_id(elems)

        if elems.startswith("is_a:"):
            parent = match_id(elems)
            edges.append((parent, uid, 'is_a'))

        if edge_type(elems):
            # Get parents stored in edge type
            eparent = match_id(elems)
            eitem = edge_type(elems)
            edges.append((eparent, uid, eitem))

        if "name:" in elems:
            name = val[1].strip().lower()

        if elems.startswith("def:"):
            definition = ":".join(val[1:]).strip()

        if uid:
            terms[uid] = [name, definition]
            names[name] = uid
            update_nodes(nodes=nodes, back_prop=back_prop, edges=edges)

    return terms, nodes, names, back_prop


def build_database():
    """
    Build the ontology database.
    """
    print(f"*** building database from: {ONTOLOGY_FILE}")

    # Check the file.
    if not os.path.isfile(ONTOLOGY_FILE):
        logger.info(f"No ontology file found, downloading...")
        download_terms()

    # Parse the terms from file
    terms, nodes, names, back_prop = parse_term(fname=ONTOLOGY_FILE)

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


def walk_tree(nodes, start, etype=None, depth=0, seen=None, collect=None):

    collect = [] if collect is None else collect
    collect.append((start, depth, etype))

    parents = nodes.get(start, [])

    seen = set()

    for par, etype in parents:

        if etype == "is_a" and par not in seen:
            walk_tree(nodes=nodes, start=par, etype=etype, seen=seen, depth=depth+1, collect=collect)

        seen.update([par])
        depth = 0


def printer(uids, terms):

    seen = set()
    for uid, etype in uids:

        if uid not in seen:
            name, definition = terms[uid]
            suffix = '' if etype == 'is_a' else f'({etype})'
            print(f"- {name} {suffix}")
        seen.update([uid])


def print_node(nodeid, terms, pad, add_highlight=False):
    vals = terms.get(nodeid, [None, None])
    name, definition = vals
    highlight = '**' if add_highlight else ''
    print(f"{pad}{nodeid}{INDENT}{highlight} {name} {highlight}")
    return


def print_tree(terms, tree=None, start=None):

    tree = [] if tree is None else tree
    tree = reversed(tree)

    # Print the definition and all children here.
    prev_depth = None
    count = 0
    for idx, objs in enumerate(tree):
        uid, depth, etype = objs

        vals = terms.get(uid, [None, None])

        if prev_depth and prev_depth == 1 and depth != 0:
            pad = INDENT * count
            count = 0
            print_node(start, terms=terms, pad=pad, add_highlight=True)
        if vals:
            add_highlight = uid.strip() == start.strip()
            pad = INDENT * count
            count += 1
            print_node(uid, terms=terms, pad=pad, add_highlight=add_highlight)

        prev_depth = depth


def show_lineage(start, terms, back_prop):

    collect = []
    walk_tree(nodes=back_prop, start=start, collect=collect)

    print_tree(tree=collect, terms=terms, start=start)

    return


def search(query, terms, prefix=""):

    # Search for names containing query.

    for uid, vals in terms.items():
        name, definition = vals

        if prefix and not uid.startswith(prefix):
            continue

        # Print all terms containing this name.
        if (query.lower() in name) or (query in uid):
            print(f"{uid}{INDENT}{name}")


def perform_query(query, terms, nodes, names, back_prop, prefix="", lineage=False):
    """
    Search database based on a name or the ontological id.
    """

    # The query is an exact match, print info.
    if names.get(query) or terms.get(query):

        # Get the GO or SO id
        uid = names.get(query) or query

        if lineage:
            show_lineage(start=uid, terms=terms, back_prop=back_prop)
            return

        # Filter for SO: or GO: ids
        if prefix and not uid.startswith(prefix):
            return

        # Get the parents.
        parents = back_prop[uid]

        # Fetch the name and definition
        name, definition = terms[uid]

        # Keep only text within quotes (if these exist).
        patt = r'"([^"]*)"'
        m = re.match(patt, definition)
        definition = m.groups()[0] if m else definition

        definition = "\n".join(wrap(definition, 80))

        print(f"\n## {name} ({uid})\n\n{definition}\n")

        print(f'Parents:')
        printer(uids=parents, terms=terms)
        children = nodes.get(uid, [])
        if children:
            print("\nChildren:")
            # Print children
            children = nodes.get(uid, [])
            printer(uids=children, terms=terms)
        print ("")

        return

    # Search for names containing query.
    search(query=query, terms=terms, prefix=prefix)


def print_stats(terms):

    gos = [k for k in terms.keys() if k.startswith('GO')]
    sos = [k for k in terms.keys() if k.startswith('SO')]
    ngos, nsos = len(gos), len(sos)
    total = ngos + nsos
    print(f"OntologyDB: total={total:,d} gene={ngos:,d} sequence={nsos:,d}")

    return


@plac.pos('query', "Search database by ontological name or GO/SO ids.")
@plac.flg('build', "build a database of all gene and sequence ontology terms. ")
@plac.flg('preload', "loads entire database in memory")
@plac.flg('verbose', "verbose mode, prints more messages")
@plac.flg('lineage', "show the ontological lineage")
@plac.flg('download', "download newest ontological data")
@plac.flg('so', "Filter query for sequence ontology terms.")
@plac.flg('go', "Filter query for gene ontology terms.")
def run(query="", build=False, download=False, preload=False, so=False, go=False,
        lineage=False, verbose=False):

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if download:
        download_terms()

    if build:
        build_database()

    terms, nodes, names, back_prop = get_data(preload=preload)

    query = query.strip()

    prefix = SO_ID if so else ''
    prefix = GO_ID if go else prefix

    if query:
        perform_query(query=query,
                      lineage=lineage,
                      terms=terms,
                      prefix=prefix,
                      nodes=nodes,
                      back_prop=back_prop,
                      names=names)
    else:
        print_stats(terms=terms)

