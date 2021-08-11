import json, shutil, os, tarfile, re
import requests
from itertools import islice
from textwrap import wrap


from biorun.libs import placlib as plac
from biorun.libs.sqlitedict import SqliteDict
from biorun import utils, const


JSON_DB = "ontology.json"
SQLITE_DB = "ontology.sqlite"


GO_URL = "http://purl.obolibrary.org/obo/go.obo"
SO_URL = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so-simple.obo"
GO_FILE = "go.obo"
SO_FILE = "so.obo"


# Delimiter in ontology files, separating one term from another.
DELIM = '[Term]'

# GO or SO id patterns.
ID_PATT = r"[G|S]O:\d+"

# Edge type pattern
EDGE_TYPE_PATT = r"(relationship:)\s+(\w+_\w+)"

GO_ID = 'GO'
SO_ID = 'SO'

join = os.path.join

# Create the full paths
GO_FILE = join(utils.DATADIR, GO_FILE)
SO_FILE = join(utils.DATADIR, SO_FILE)
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


def download_prebuilt():
    utils.download_from_bucket(bucket_name=const.BUCKET_NAME, file_name=JSON_DB, cache=True)
    utils.download_from_bucket(bucket_name=const.BUCKET_NAME, file_name=SQLITE_DB, cache=True)


def download_terms():

    utils.download(url=GO_URL, dest_name=GO_FILE)
    utils.download(url=SO_URL, dest_name=SO_FILE)

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
        terms = utils.open_db(TERM, fname=SQLITE_DB)
        nodes = utils.open_db(GRAPH, fname=SQLITE_DB)
        names = utils.open_db(NAMES, fname=SQLITE_DB)
        back = utils.open_db(CHILDREN, fname=SQLITE_DB)

    return terms, nodes, names, back


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
    etype = match.group(2) if match else None

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
    print(f"### parsing: {fname}")
    for elems in stream:

        if stop_term(elems):
            # Reset objs for next term.

            if uid:
                terms[uid] = [name, definition]
                names[name] = uid
                update_nodes(nodes=nodes, back_prop=back_prop, edges=edges)
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

    return terms, nodes, names, back_prop


def build_database(fname, flg='w'):
    """
    Build the ontology database.
    """
    print(f"### building database from: {fname}")

    # Check the file.
    if not os.path.isfile(fname):
        logger.info(f"No ontology file found, downloading...")
        download_terms()

    # Parse the terms from file
    terms, nodes, names, back_prop = parse_term(fname=fname)

    save = lambda table, vals:  utils.save_table(table, vals, fname=SQLITE_DB, flg=flg)

    # Save terms into the database
    save(TERM, terms)

    # Save graph into the database
    save(GRAPH, nodes)

    # Save names into the database
    save(NAMES, names)

    # Save child-parent into the database
    save(CHILDREN, back_prop)

    return terms, nodes, names, back_prop


def build_db():
    """
    Wrapper to build both GO and SO.
    """
    terms, nodes, names, back_prop = build_database(GO_FILE)
    # Add to the JSON model, instead of rewrtiting it.
    soterms, sonodes, sonames, soback_prop = build_database(SO_FILE, flg='c')

    terms.update(soterms)
    nodes.update(sonodes)
    names.update(sonames)
    back_prop.update(soback_prop)

    print("### saving the JSON model")
    store = dict(TERMS=terms, GRAPH=nodes, NAMES=names, CHILDREN=back_prop)
    fp = open(JSON_DB, "wt")
    json.dump(store, fp, indent=4)
    fp.close()

    return


def walk_tree(nodes, start, etype=None,
              visited=None,
              collect=None,
              one_pass=False):

    collect = [] if collect is None else collect
    visited = set() if visited is None else visited
    collect.append(start)

    parents = nodes.get(start, [])

    for node, et in parents:
        if node in visited:
            continue

        if (et == etype) or (etype is None):
            visited.update([node])
            walk_tree(nodes=nodes,
                      start=node,
                      etype=etype,
                      visited=visited,
                      collect=collect,
                      one_pass=one_pass)
        if one_pass:
            break


def printer(uids, terms, prefix=""):

    seen = set()
    for uid, etype in uids:

        if uid not in seen:
            name, definition = terms[uid]
            suffix = '' if etype == 'is_a' else f'({etype})'
            print(f"{prefix}- {name} {suffix}")
        seen.update([uid])


def wrap_text(text, pre=''):
    # Keep only text within quotes (if these exist).
    patt = r'"([^"]*)"'
    m = re.match(patt, text)
    text = m.groups()[0] if m else text
    text = f"\n{pre}".join(wrap(text, 80))
    return text


def formatted_printer(name, uid, definition,  nodes, terms, parents=[], prefix=""):
    definition = wrap_text(definition, pre=prefix)
    print(f"\n{prefix}## {name} ({uid})\n\n{prefix}{definition}\n")
    if parents:
        print(f'{prefix}Parents:')
        printer(uids=parents, terms=terms, prefix=prefix)
    children = nodes.get(uid, [])

    if children:
        print(f"\n{prefix}Children:")
        # Print children
        printer(uids=children, terms=terms, prefix=prefix)
    print(" ")
    return


def print_node(uid, name, nodes=None, terms=None, define='', pad='', is_leaf=False):

    if is_leaf:
        formatted_printer(nodes=nodes, terms=terms,
                          name=name,
                          uid=uid, definition=define,
                          parents=[], prefix=pad)
    else:
        print(f"{pad}{uid}{INDENT}{name}")

    return


def show_lineage(start, terms, nodes, back_prop):

    collect = []

    walk_tree(nodes=back_prop, start=start, collect=collect, one_pass=True, etype='is_a')

    tree = reversed(collect)

    # Print the definition and all children here.
    for depth, uid in enumerate(tree):
        vals = terms.get(uid)
        if vals:
            names, define = vals
            is_leaf = uid.strip() == start.strip()
            pad = INDENT * depth
            print_node(uid=uid, name=names,
                       define=define,
                       pad=pad, is_leaf=is_leaf,
                       nodes=nodes, terms=terms)

    parents = back_prop.get(start, [])
    if len(parents) > 1:
        print("### More than on path detected, use -P to view all relationships.\n")

    return


def search(query, terms, prefix=""):

    # Search for names containing query.

    for uid, vals in terms.items():
        name, definition = vals

        if prefix and not uid.startswith(prefix):
            continue

        # Print all terms containing this name.
        if (query.lower() in name) or (query in uid):
            print_node(uid=uid, name=name)


def perform_query(query, terms, nodes, names, back_prop, prefix="", lineage=False):
    """
    Search database based on a name or the ontological id.
    """

    # The query is an exact match, print info.
    if names.get(query) or terms.get(query):

        # Get the GO or SO id
        uid = names.get(query) or query

        if lineage:
            show_lineage(start=uid, terms=terms, back_prop=back_prop, nodes=nodes)
            return

        # Filter for SO: or GO: ids
        if prefix and not uid.startswith(prefix):
            return

        # Get the parents.
        parents = back_prop.get(uid,[])

        # Fetch the name and definition
        name, definition = terms[uid]

        formatted_printer(name=name, uid=uid,
                          definition=definition,
                          parents=parents, nodes=nodes,
                          terms=terms)

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


def plot_term(query, names, terms, nodes, back_prop, outname=''):

    try:
        import pygraphviz as pgv
    except ImportError as exc:
        utils.error(exc, stop=False)
        utils.error("Try: conda install pygraphviz")

    if names.get(query) or terms.get(query):
        uid = names.get(query) or query
    else:
        return

    def frmt(uid, name):
        out = fr"{uid}\n{name}"
        return out

    collect = []

    walk_tree(nodes=back_prop, start=uid, collect=collect)

    grph = pgv.AGraph(directed=True)

    # Iterate through nodes and build tree.
    for item in collect:
        # Get all children for this node present in tree
        children = nodes.get(item, [])
        chls = [c for c in children if c[0] in collect]
        name, define = terms.get(item)

        # Pair item with each child as an edge.
        for child in chls:
            chl, etype = child
            cname, cdefine = terms.get(chl)

            style = "solid" if etype == "is_a"else "dashed"

            # Format the edge to include both id and name.
            grph.add_edge(frmt(item, name), frmt(chl, cname), label=f" {etype}", style=style)

        # Add a leaf node
        if not children:
            grph.add_node(frmt(item, name))

    grph.edge_attr.update(shape="normal", color='black', dir="back", penwidth=2)
    grph.node_attr.update(shape="box", style="rounded,filled", fillcolor="beige")

    # Highlight the query term.
    name, define = terms.get(uid)
    try:
        node = grph.get_node(frmt(uid, name))
        node.attr.update(fillcolor="plum")
    except Exception as exc:
        logger.error(exc)
        pass

    # Construct file name and write to pdf.
    grph.layout(prog='dot')

    print(f"### Writing plot to {outname}")
    # Write DOT string to file and plot to .pdf
    if outname.endswith('.dot'):
        open(outname, 'w').write(grph.to_string())
    else:
        grph.draw(outname)

    return


@plac.pos('query', "Search database by ontological name or GO/SO ids.")
@plac.flg('build', "build a database of all gene and sequence ontology terms. ")
@plac.flg('preload', "loads entire database in memory", abbrev="P")
@plac.flg('verbose', "verbose mode, prints more messages")
@plac.flg('lineage', "show the ontological lineage")
@plac.flg('download', "download prebuilt database")
@plac.flg('so', "Filter query for sequence ontology terms.")
@plac.flg('go', "Filter query for gene ontology terms.")
@plac.flg('update', "Update latest terms from remote hosts and build with those.")
@plac.opt('plot', "Plot the network graph of the given GO term into the given file name.", abbrev="p")
@plac.flg('define', "search ontology for definitions", abbrev='x')
def run(build=False, download=False, preload=False, so=False, go=False,
        lineage=False, update=False, plot='', define=False, verbose=False, *query):

    # Join up all words.
    query = " ".join(query)

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    if download:
        download_prebuilt()

    if update:
        download_terms()

    if build:
        build_db()

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

    if plot:
        plot_term(query=query, names=names, terms=terms, nodes=nodes, back_prop=back_prop, outname=plot)
