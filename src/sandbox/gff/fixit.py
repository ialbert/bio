"""
Adds version number sequence id of a GFF file.
"""
import sys
ACC=sys.argv[1]
VER=sys.argv[2]
for line in sys.stdin:
    line = line.strip()
    elems = line.split()
    if elems and elems[0] == ACC:
        elems[0] = f'{elems[0]}.{VER}'

    if line.startswith("#"):
        print (line)
    else:
        print("\t".join(elems))