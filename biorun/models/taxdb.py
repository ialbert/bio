import codecs
import csv
import json
import os
import re
import sys
import tarfile
from itertools import islice

from biorun import fetch, const, ncbi
from biorun import utils
from biorun.libs import placlib as plac
from biorun.libs.sqlitedict import SqliteDict
from biorun.models import jsonrec

JSON_DB_NAME = "taxdb.json"
SQLITE_DB_NAME = "taxdb.sqlite"
TAXDB_URL = "http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAXDB_NAME = "taxdump.tar.gz"

join = os.path.join

# Create the full paths
TAXDB_NAME = join(utils.DATADIR, TAXDB_NAME)
SQLITE_DB = join(utils.DATADIR, SQLITE_DB_NAME)
JSON_DB = join(utils.DATADIR, JSON_DB_NAME)

# Create the thing here.

GRAPH, BACK, NAMES, SYNON, LATIN = "GRAPH", "BACKLINKS", "NAMES", "SYNONYMS", "LATIN"

# Indentation character
INDENT = '  '

# Fields separator

SEP = ', '

# Used during debugging only to speed up database builds.
# Keep it at None
LIMIT = None

CHUNK = 25000

logger = utils.logger


def download_prebuilt():
    """
    Download prebuild databases.
    """
    utils.download_from_bucket(bucket_name=const.BUCKET_NAME, file_name=SQLITE_DB_NAME, cache=True)
    utils.download_from_bucket(bucket_name=const.BUCKET_NAME, file_name=JSON_DB_NAME, cache=True)

    # Download the taxonomy file.
    update_taxdump()


def update_taxdump(url=TAXDB_URL, dest_name=TAXDB_NAME):
    """
    Downloads taxdump file.
    """
    utils.download(url=url, dest_name=dest_name)


def get_stream(tar, name, limit=None):
    mem = tar.getmember(name)
    stream = tar.extractfile(mem)
    stream = map(lambda x: x.decode("ascii"), stream)
    stream = islice(stream, limit)
    return stream


def search_names(word, fname=TAXDB_NAME, name="names.dmp", limit=None):
    """
    Processes the names.dmp component of the taxdump.
    """

    # Needs a taxdump to work.
    if not os.path.isfile(fname):
        utils.error("taxdump file not found (download and build it first)")

    # The taxdump file.
    tar = tarfile.open(fname, "r:gz")
    stream = get_stream(tar=tar, name=name, limit=limit)
    stream = csv.reader(stream, delimiter="\t")

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


def parse_names(fname, name="names.dmp", limit=None):
    """
    Parses the names.dmp component of the taxdump.
    """

    # The taxdump file.
    tar = tarfile.open(fname, "r:gz")

    stream = get_stream(tar=tar, name=name, limit=limit)
    stream = csv.reader(stream, delimiter="\t")

    # Various lookup tables.
    tax2data, name2tax = {}, {}

    print(f"*** parsing: {name}")

    for index, elems in enumerate(stream):
        taxid, name, label = elems[0], elems[2], elems[6]

        # Assembly count (TODO)
        count = 0

        # Scientific name fields.
        if label == 'scientific name':
            # Name, rank, common name, parent, assembly count
            tax2data[taxid] = [name, "", "", "", count]
        elif label == 'genbank common name':
            name2tax[name] = taxid

    # Backfill common names when these are different.
    for taxid, name in name2tax.items():
        if taxid in tax2data:
            tax2data[taxid][2] = name

    return tax2data


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

    # Get the assembly.
    # _, _, taxon_acc = ncbi.parse_summary()

    # Parse the names
    name_dict = parse_names(fname, limit=limit)

    # Parse the nodes.
    node_dict, back_dict = parse_nodes(fname, name_dict=name_dict, limit=limit)

    def save_table(name, obj):
        utils.save_table(name=name, obj=obj, fname=SQLITE_DB)

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


def print_metadata(terms):
    def formatter(row):
        print("\t".join(row))

    for term in terms:
        lines = get_metadata(term)
        lines = filter(lambda x: x.split(), lines)
        old_header = next(lines)
        new_header = "host species accession date location isolate".split()

        print("\t".join(new_header))
        for line in lines:
            print(line)


def get_metadata(taxid, limit=None):
    """
    Returns all accessions
    """
    import requests

    # The dataset accession point.
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v1alpha/virus/taxon/{taxid}/genome/table"

    params = {
        'format': 'tsv',
        'refseq_only': "false",
        'complete_only': 'true',
        'table_fields': [
            'host_tax_id', 'species_tax_id',
            'nucleotide_accession',
            'collection_date', 'geo_location', 'isolate_name',
        ]
    }

    conn = requests.get(url, stream=True, params=params)
    lines = conn.iter_lines()
    lines = islice(lines, limit)

    if conn.status_code != 200:
        msg = f"HTTP status code: {conn.status_code}"
        utils.error(msg)

    lines = map(decode, lines)

    return lines


def get_values(node, names):
    # Avoiding code duplication everywhere
    sname, rank, cname, parent, acount = names.get(node, ("MISSING", "NO RANK", "", "", 0))
    return sname, rank, cname, parent, acount


def print_assemblies(taxid, assembly):
    """
    Print assemblies
    """

    assemblies = assembly.get(taxid, [])

    for acc in assemblies:
        print(f'{taxid}{INDENT}{acc}')


def node_formatter(node, names, depth):
    """
    Creates a long form representation of a node.
    """
    indent = INDENT * depth
    sname, rank, name, parent, count = get_values(node, names)

    # Get any full genome assemblies this node may have.

    # Decide what to do with common names.
    if name and name != sname:
        data = [rank, node, sname, name]
    else:
        data = [rank, node, sname]

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


def get_data(preload=False, acc=False):
    """
    Returns the graph structure for the database.
    """
    if preload:
        if not os.path.isfile(JSON_DB):
            utils.error(f"taxonomy file not found (you must build it first): {JSON_DB}")
        store = json.load(open(JSON_DB))
        names = store[NAMES]
        graph = store[GRAPH]
        latin = store[LATIN]
    else:
        names = open_db(NAMES)
        graph = open_db(GRAPH)
        latin = open_db(LATIN)

    if acc:
        _, taxon_acc, _ = ncbi.get_data()
    else:
        taxon_acc = {}

    return names, graph, taxon_acc, latin


def print_stats(names, graph):
    node_size, graph_size = len(names), len(graph)
    print(f"TaxDB: nodes={node_size:,d} parents={graph_size:,d}")


def search_taxa(word, preload=False):
    names, graph, assembly, latin = get_data(preload=preload)

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


def filter_file(stream, keep, remove, graph, colidx=0):
    """
    Filters a file to retain only the rows where a taxid is ina subtree.
    """
    # Collects all children of the taxids.
    keep_dict, remove_dict = {}, {}

    # Collect the matching nodes.
    for term in keep.split(","):
        dfs_visitor(graph=graph, node=term, visited=keep_dict)

    # Collect the matching nodes.
    for term in remove.split(","):
        dfs_visitor(graph=graph, node=term, visited=remove_dict)

    # Read the stream.
    reader = csv.reader(stream, delimiter="\t")

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
    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerows(reader)


def parse_taxids(json):
    """
    Attempts to parse taxids from a json data
    """
    # Parses the taxids
    doubles = [jsonrec.find_taxid(rec) for rec in json] if json else [[]]
    # Flatten the list
    taxids = [elem for sublist in doubles for elem in sublist]
    return taxids


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


def parse_lines(stream, field=1, sep="\t"):
    colidx = field - 1
    stream = filter(lambda x: len(x.split(sep)) >= colidx, stream)
    stream = filter(lambda x: not x.startswith('#'), stream)
    return stream


@plac.pos("terms", "taxids or search queries")
@plac.flg('update', "updates and builds a local database")
@plac.flg('preload', "loads entire database in memory")
@plac.flg('list_', "lists database content", abbrev='l')
@plac.opt('scinames', "scientific or common names in each line. ", abbrev="S")
@plac.flg('children', "include children when returning when parsing latin names", abbrev='C')
@plac.flg('lineage', "show the lineage for a taxon term", abbrev="L")
@plac.opt('indent', "the indentation depth (set to zero for flat)")
@plac.opt('sep', "separator (default is ', ')", abbrev='s')
@plac.flg('metadata', "downloads metadata for the taxon", abbrev='m')
@plac.flg('download', "downloads the database from the remote site", abbrev='G')
@plac.opt('depth', "how deep to visit a clade ", abbrev='d', type=int)
@plac.opt('keep', "clade to keep", abbrev='K')
@plac.opt('remove', "clade to remove", abbrev='R')
@plac.opt('field', "which column to read when filtering")
@plac.flg('verbose', "verbose mode, prints more messages")
@plac.flg('accessions', "Print the accessions number for each ")
def run(lineage=False, update=False, download=False, accessions=False, keep='', remove='', field=1,
        scinames='', children=False, list_=False, depth=0, metadata=False, preload=False, indent=2, sep='',
        verbose=False, *terms):
    global SEP, INDENT, LIMIT

    # Input connected to a stream
    if not sys.stdin.isatty():
        terms = sys.stdin.readlines()

    # Indentation level
    INDENT = ' ' * indent

    # Separator string.
    SEP = decode(sep) if sep else ", "

    # Set the verbosity
    utils.set_verbosity(logger, level=int(verbose))

    # Access the database.
    names, graph, assembly, latin = get_data(preload=preload, acc=accessions)

    # Download prebuilt database.
    if download:
        download_prebuilt()

    # Updates the taxdump and builds a new taxonomy file.
    if update:
        # update_taxdump()
        build_database(limit=LIMIT)

    # List the content of a database.
    if list_:
        print_database(names=names, graph=graph)
        sys.exit()

    # Obtain metadata for the taxon
    if metadata:
        print_metadata(terms)
        return

    if scinames:
        search_file(scinames, names=names, latin=latin, graph=graph, include=children)
        sys.exit()

    # Filters a file by a column.
    if keep or remove:
        filter_file(stream=terms, keep=keep, remove=remove, graph=graph, colidx=field - 1)
        sys.exit()


    # Input may come from a file or command line.
    colidx = field - 1
    terms = parse_lines(terms)
    terms = map(lambda x: x.split("\t")[colidx].strip(), terms)
    terms = filter(None, terms)
    terms = list(terms)

    # No valid terms found. Print database stats.
    if not terms:
        print_stats(names=names, graph=graph)
        sys.exit()

    # These are the terms looked up in the database.
    words = []

    # Some terms may be valid data names.
    for term in terms:
        term = term.strip()
        # Attempts to interpret the word as an existing dataset.
        json = fetch.get_json(term)

        # Extend the search temrs.
        taxids = parse_taxids(json) if json else [term]

        # Add to the terms.
        words.extend(taxids)

    # Produce lineages
    if lineage:
        for term in words:
            print_lineage(term, names=names)
        sys.exit()

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
