import codecs
import csv
import json
import os
import re
import sys
import tarfile
from itertools import islice

from biorun import convert
from biorun import utils
from biorun.libs import placlib as plac
from biorun.libs.sqlitedict import SqliteDict

# Database filenames.
SQLITE_FILE = "taxonomy.sqlite"
SQLITE_PATH = utils.cache_path(SQLITE_FILE)

JSON_FILE = "taxonomy.json"
JSON_PATH = utils.cache_path(JSON_FILE)

# NCBI taxonomy files.
TAXDB_FILE = "taxdump.tar.gz"
TAXDB_PATH = utils.cache_path(TAXDB_FILE)
TAXDB_URL = "http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Keys into the database
GRAPH, BACK, TAXID = "GRAPH", "BACK", "TAXIDS"

# Indentation character
INDENT = '  '

# Fields separator
SEP = ', '

# Used during debugging only to speed up database builds.
# Keep it at None
LIMIT = None

CHUNK = 25000

logger = utils.logger

def build_database(limit=None):
    """
    Downloads taxdump file.
    """

    # The location of the taxdump file.
    path = TAXDB_PATH

    # Download the latest taxdump file.
    utils.download(url=TAXDB_URL, fname=path)

    print(f"\n# taxdump file: {path}")

    # Parse the names
    tax2data = parse_names(path, limit=limit)

    # Parse the nodes and backpropagation.
    graph = parse_nodes(path, tax2data=tax2data, limit=limit)

    # A shortcut to the function.
    def save_table(name, obj):
        utils.save_table(name=name, obj=obj, fname=SQLITE_PATH)

    # Save the taxid definitions.
    save_table(TAXID, tax2data)

    # Save the graph.
    save_table(GRAPH, graph)

    print("# saving the JSON model")
    json_path = os.path.join(utils.DATADIR, JSON_PATH)

    # Save the JSON file as well.
    store = {TAXID: tax2data, GRAPH: graph}
    fp = open(json_path, 'wt')
    json.dump(store, fp, indent=4)
    fp.close()


def open_tarfile(archive, filename, limit=None, delimiter="\t"):
    """
    Returns the content of named file in a tarred archive.
    """
    tar = tarfile.open(archive, "r:gz")
    mem = tar.getmember(filename)
    stream = tar.extractfile(mem)
    stream = map(lambda x: x.decode("ascii"), stream)
    stream = islice(stream, limit)
    stream = csv.reader(stream, delimiter=delimiter)
    return stream


def search_names(word, archive=TAXDB_PATH, name="names.dmp", limit=None):
    """
    Processes the names.dmp component of the taxdump.
    """

    # Needs a taxdump to work.
    if not os.path.isfile(archive):
        utils.error("taxdump file not found (download and build it first)")

    # Open stream into the tarfile.
    stream = open_tarfile(archive=archive, filename=name, limit=limit)

    # The pattern may be regular expression.
    patt = re.compile(word, re.IGNORECASE)

    # Labels that will be searched.
    valid = {'scientific name', 'equivalent name', 'genbank common name'}

    def select(row):
        taxid, name, label = row[0], row[2], row[6]
        return label in valid and patt.search(name)

    # Apply the selector.
    stream = filter(select, stream)
    for elems in stream:
        taxid, name, label = elems[0], elems[2], elems[6]
        yield taxid, name


def parse_names(archive, filename="names.dmp", limit=None):
    """
    Parses the names.dmp component of the taxdump.
    """

    # Parse the tarfile.
    stream = open_tarfile(archive=archive, filename=filename, limit=limit)

    # Lookup tables.
    tax2data, name2tax = {}, {}

    print(f"# processing: {filename}")

    # Process the nodes.dmp file.
    for row in stream:

        # The names.dmp file structure.
        taxid, name, label = row[0], row[2], row[6]

        # Assembly count if exists (TODO)
        count = 0

        # Populate only for scientific names.
        if label == 'scientific name':
            # 5 columns: sciname, rank, common name, parent, assembly count
            # Some information will be only known later, from the nodes.dmp
            tax2data[taxid] = [name, "", "", "", count]
        elif label == 'genbank common name':
            name2tax[name] = taxid

    # Fill common names when it exists.
    for taxid, name in name2tax.items():
        if taxid in tax2data:
            tax2data[taxid][2] = name

    return tax2data


def parse_nodes(archive, tax2data, filename="nodes.dmp", limit=None):
    """
    Parses the names.dmp component of the taxdump.
    """

    # Parse the NCBI taxump.
    stream = open_tarfile(archive=archive, filename=filename, limit=limit)

    # Data structures to fill.
    graph = {}

    print("# processing: nodes.dmp")

    # Process the nodes.dmp file.
    for row in stream:

        # nodes.dmp file format.
        child, parent, rank = row[0], row[2], row[4]

        # Connect parent to all children
        graph.setdefault(parent, []).append(child)

        # Mutates existing datastructure with rank and parent info.
        if child in tax2data:
            tax2data[child][1] = rank
            tax2data[child][3] = parent

    return graph


queue = list()


def get_values(node, names):
    sciname, rank, name, parent, count = names.get(node, ("MISSING", "NO RANK", "", "", 0))
    return sciname, rank, name, parent, count


def node_formatter(node, names, depth):
    """
    Creates a long form representation of a node.
    """
    indent = INDENT * depth
    sciname, rank, name, parent, count = get_values(node, names)

    # Write common name if exists and different from sciname
    if name and name != sciname:
        data = [rank, node, sciname, name]
    else:
        data = [rank, node, sciname]

    text = indent + SEP.join(data)

    return text


def backprop(node, names, collect=[]):
    """
    Collects nodes when propagating backwards.
    """
    if node in names:
        sciname, rank, name, parent, count = names[node]
        if parent and parent != node:
            collect.append(parent)
            backprop(parent, names, collect)


def print_lineage(taxid, names, flat=0):
    """
    Prints the lineage for a taxid.
    """

    # Must be a valid taxid.
    if taxid not in names:
        msg = f"Invalid taxid: {taxid}"
        utils.error(msg)

    # Will back propagate to parents.
    collect = [taxid]
    backprop(taxid, names, collect=collect)

    # Going back to superkingdom only.
    collect = collect[:-1]

    # Start at the parent.
    collect = reversed(collect)

    # Format each node.
    for step, node in enumerate(collect):
        text = node_formatter(node, names=names, depth=step)
        print(text)


def get_data(preload=False, strict=True):
    """
    Returns the graph structure for the database.
    """
    if preload:
        if not os.path.isfile(JSON_PATH):
            utils.error(f"taxonomy file not found (you must build it first): {JSON_PATH}")
        store = json.load(open(JSON_PATH))
        names = store[TAXID]
        graph = store[GRAPH]
    else:
        names = utils.open_db(TAXID, fname=SQLITE_PATH, strict=strict)
        graph = utils.open_db(GRAPH, fname=SQLITE_PATH, strict=strict)

    return names, graph


def print_stats(names, graph):
    node_size, graph_size = len(names), len(graph)
    print(f"Content: {node_size:,d} taxon names; {graph_size:,d} parent ranks")


def search_taxa(word, preload=False):
    names, graph = get_data(preload=preload)

    word = decode(word)

    print(f"# Searching taxonomy for: {word}")
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


def print_term(taxid, graph, names, maxdepth=0):
    """
    Prints a term when visited via DFS.
    """

    def formatter(node, depth, **kwds):
        text = node_formatter(node, names=names, depth=depth)
        print(text)

    dfs_visitor(graph, taxid, visited={}, func=formatter, maxdepth=maxdepth)


def donothing(*args, **kwds):
    """
    Placeholder to perform no action.
    """
    pass


def valid_int(text):
    try:
        int(text)
        return True
    except ValueError as exc:
        return False


def dfs_visitor(graph, node, visited, depth=0, func=donothing, maxdepth=0):
    """
    Performs depth-first search and collects output into the visited dictionary keyed by depth.
    Calls func at every visit.
    """
    if node not in visited:
        visited[node] = depth
        func(node=node, depth=depth, visited=visited)
        for nbr in graph.get(node, []):
            nextdepth = depth + 1
            if maxdepth and nextdepth >= maxdepth:
                continue
            dfs_visitor(graph=graph, node=nbr, depth=nextdepth, visited=visited, func=func, maxdepth=maxdepth)


def filter_file(stream, terms, keep, remove, graph, colidx=0, sep="\t"):
    """
    Filters a file to retain only the rows where a taxid is ina subtree.
    """
    if not stream:
        if len(terms) == 0:
            msg = f"filtering needs an input stream or a filename"
            utils.error(msg)
        stream = open(terms[0])

    # Collects all children of the taxids.
    keep_dict, remove_dict = {}, {}

    # Taxids to keep
    keeper = keep.split(",")
    # Fill the keeper dictionary.
    for term in keeper:
        dfs_visitor(graph=graph, node=term, visited=keep_dict)

    # Fill the remover dictionary.
    remover = remove.split(",")
    for term in remover:
        dfs_visitor(graph=graph, node=term, visited=remove_dict)

    # Read the stream.
    reader = csv.reader(stream, delimiter=sep)

    # Selection condition.
    def keep_func(row):
        taxid = row[colidx]
        return taxid in keep_dict

    def remove_func(row):
        taxid = row[colidx]
        return taxid not in remove_dict

    # What to keep.
    if keep:
        reader = filter(keep_func, reader)

    # What to remove.
    if remove:
        reader = filter(remove_func, reader)

    # Generate the output.
    writer = csv.writer(sys.stdout, delimiter=sep)
    writer.writerows(reader)



def decode(text):
    """
    Recognize string encodings: \t etc
    """
    return codecs.decode(text, 'unicode_escape')


def isnum(x):
    try:
        int(x)
        return True
    except ValueError as exc:
        return False


def parse_stream(stream, field=1, delim="\t"):
    """
    Parses a stream for a column.
    """
    # One based coordinate system.
    colidx = field - 1

    # Sanity check.
    assert colidx >= 0

    # Remove comments
    stream = filter(lambda x: not x.startswith('#'), stream)

    # Remove empty lines
    stream = filter(lambda x: x.strip(), stream)

    # Create a reader.
    reader = csv.reader(stream, delimiter=delim)

    # Keep only rows that have data for the column
    reader = filter(lambda row: len(row) >= colidx, reader)

    # Exract the correct column
    reader = (row[colidx] for row in reader)

    # Keep only nonzero entries.
    reader = filter(lambda x: x.strip(), reader)

    return list(reader)


@plac.pos("terms", "taxids or search queries")
@plac.flg('build', "updates and builds a local database")
@plac.flg('preload', "loads entire database in memory")
@plac.flg('list_', "lists database content", abbrev='l')
# @plac.opt('scinames', "scientific or common names in each line. ", abbrev="S")
# @plac.flg('children', "include children when returning when parsing latin names", abbrev='C')
@plac.flg('lineage', "show the lineage for a taxon term", abbrev="L")
@plac.opt('indent', "the indentation depth (set to zero for flat)")
@plac.opt('sep', "separator (default is ', ')", abbrev='s')
@plac.opt('depth', "how deep to visit a clade ", abbrev='d', type=int)
@plac.opt('keep', "clade to keep", abbrev='K')
@plac.opt('remove', "clade to remove", abbrev='R')
@plac.opt('field', "which column to read when filtering")
@plac.flg('verbose', "verbose mode, prints more messages")
@plac.flg('accessions', "Print the accessions number for each ")
def run(lineage=False, build=False, accessions=False,
        keep='', remove='', field=1,
        list_=False, depth=0, preload=False, indent=2, sep=',',
        verbose=False, *terms):
    global SEP, INDENT, LIMIT

    terms = list(terms)

    # Do not consume stdin when filtering input
    filter_mode = keep or remove

    # Input may come as a stream.
    if not sys.stdin.isatty() and not filter_mode:
        ids = parse_stream(sys.stdin, field=field, delim=sep)
        # Dictionaries maintain the order.
        store = dict((map(lambda x: (x, x), ids)))
        terms.extend(store.keys())

    # Indentation level
    INDENT = ' ' * indent

    # Separator string.
    SEP = decode(sep) if sep else ", "

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Downloads a new taxdump and builds a new taxonomy database.
    if build:
        build_database(limit=LIMIT)
        return

    # Get the content of the database.
    names, graph = get_data(preload=preload)

    # List the content of a database.
    if list_:
        print_database(names=names, graph=graph)
        return



    # if scinames:
    #    search_file(scinames, names=names, latin=latin, graph=graph, include=children)
    #    return

    # Filters a file by a column.
    if keep or remove:
        filter_file(stream=sys.stdin, terms=terms, keep=keep, remove=remove, graph=graph, colidx=field-1, sep=sep)
        return

    # No valid terms found. Print database stats.
    if not terms:
        print_stats(names=names, graph=graph)
        return

    # These are the terms looked up in the database.
    words = []

    # Some terms may be valid data names.
    for term in terms:
        term = term.strip()

        if os.path.isfile(term):
            # Attempts to parse taxonomical information from GenBank
            recs = convert.parse_stream(fname=term)
            recs = filter(convert.source_only(True), recs)
            for rec in recs:
                db_xref = rec.annot.get('db_xref', [])
                db_xref = filter(lambda x: x.startswith("taxon:"), db_xref)
                db_xref = map(lambda x: x.split(":")[1], db_xref)
                db_xref = list(db_xref)
                words.extend(db_xref)
        else:
            # Add to the terms.
            words.append(term)

    # Produce lineages
    if lineage:
        for term in words:
            print_lineage(term, names=names)
        return

    # Will check to mixed terms (valid taxids and search words mixed)

    # Truth vector to terms in names.
    valid = list(map(lambda x: x in names, words))
    any_valid = any(valid)
    all_valid = all(valid)

    # Mixed term condition.
    mixed_terms = any_valid and not all_valid

    # We don't allow mixed terms (produces different outputs).
    if mixed_terms:
        invalid = ", ".join(filter(lambda x: x not in names, words))
        msg = f"Unkown taxids: {invalid}"
        utils.error(msg)

    # Apply the approprate task to each term separately.
    for term in words:
        if all_valid:
            print_term(term, names=names, graph=graph, maxdepth=depth)
        else:
            search_taxa(term)


if __name__ == '__main__':
    # Bony fish: 117565
    # Betacoronavirus: 694002
    # SARS-COV2: 2697049

    # Jawless vertebrates: 1476529
    # Vertebrata: 7742

    plac.call(run)
