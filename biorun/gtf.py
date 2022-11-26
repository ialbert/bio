"""
TODO:

Generates tab delimited files that represent mappings between attributes of a GTF file.

For example:

transcript_id  gene_id

"""

USAGE = f'''

    cat genes.gtf | bio gtf


'''
from itertools import *
import csv, sys
import biorun.libs.placlib as plac

def parse(ftype, source, target, limit):

    limit = None if not limit else int(limit)

    ftype = '' if ftype == 'all' else ftype

    # Valid GFF files have 9 columns
    def valid(x):
        return len(x) == 9

    # Filter for features of a given type.
    def keeper():

        def typer(x):
            return x[2] == ftype

        def ident(x):
            return True

        func = typer if ftype else ident
        return func

    # Split and process attributes.
    def split_attrs(elems):
        attrs = elems[8].strip(' ;')
        attrs = attrs.split(";")
        attrs = [field.split() for field in attrs]
        attrs = [(field[0], field[1].strip(' "')) for field in attrs ]
        return attrs

    # Build the
    stream = sys.stdin
    stream = csv.reader(stream, delimiter="\t")
    stream = islice(stream, limit)
    stream = filter(valid, stream)
    stream = filter(keeper(), stream)
    stream = map(split_attrs, stream)

    # Extract mapping for each feature
    for index, elems in enumerate(stream):

        row = dict(elems)
        try:
            source_id = row[source]
        except KeyError as exc:
            msg1 = f'# {elems}\n'
            msg2 = f'# Error: missing {exc} key at line {index+1}\n'
            sys.stderr.write(msg1)
            sys.stderr.write(msg2)
            sys.exit(1)

        target_id = row.get(target, "Unknown")
        data = [source_id, target_id]

        yield data

def combine():
    pass

@plac.opt('type_', "GTF type filter", abbrev="T", type=str)
@plac.opt('source', "source attribute", abbrev="s", type=str)
@plac.opt('target', "target attribute", abbrev="t", type=str)
@plac.opt('limit', "how many lines to parse", abbrev="l", type=int)
@plac.flg('collate', "how many lines to parse", abbrev="l", type=int)
def run(type_='transcript', source="transcript_id", target="gene_id", limit=None, collate=False):
    """
    cat genes.gtf | bio gtf

    cat genes.gtf | bio gtf --target gene_name

    cat genes.gtf | bio gtf --type gene

    """
    stream = parse(ftype=type_, source=source, target=target, limit=limit)

    if collate:
        pass

    for row in stream:
        print("\t".join(row))

def main():
    """
    Entry point for the script.
    """
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    plac.call(run)

if __name__ == '__main__':
    main()

