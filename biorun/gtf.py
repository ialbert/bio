"""
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

    def valid(x):
        return len(x) == 9

    def keeper():
        def func(x):
            return x[2] == ftype

        return func

    def splitter(elems):
        attrs = elems[8].strip(' ;')
        attrs = attrs.split(";")
        attrs = [field.split() for field in attrs]
        attrs = [(field[0], field[1].strip(' "')) for field in attrs ]
        return attrs

    stream = sys.stdin
    stream = csv.reader(stream, delimiter="\t")
    stream = islice(stream, limit)
    stream = filter(valid, stream)
    if ftype:
        stream = filter(keeper(), stream)
    stream = map(splitter, stream)
    for index, elems in enumerate(stream):
        #print (elems)
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
        print("\t".join(data))


@plac.opt('type_', "GTF type filter", abbrev="T", type=str)
@plac.opt('source', "source attribute", abbrev="s", type=str)
@plac.opt('target', "target attribute", abbrev="t", type=str)
@plac.opt('limit', "how many lines to parse", abbrev="l", type=int)
def run(type_='transcript', source="transcript_id", target="gene_id", limit=None):
    """
    cat genes.gtf | bio gtf

    cat genes.gtf | bio gtf --target gene_name

    cat genes.gtf | bio gtf --type gene

    """
    parse(ftype=type_, source=source, target=target, limit=limit)

def main():
    """
    Entry point for the script.
    """
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    plac.call(run)

if __name__ == '__main__':
    main()

