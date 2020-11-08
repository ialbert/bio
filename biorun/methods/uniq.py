import fileinput
from collections import defaultdict

def run():
    for line in fileinput.input():
        line = line.strip()
        if line.startswith("#"):
            continue
        if line:
            pass