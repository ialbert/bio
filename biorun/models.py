import sys
from itertools import *

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Alignment:
    """
    Represents a pairwise alignment.
    """

    def __init__(self, target: SeqRecord, query: SeqRecord, score=None):

        # Aligned query and target.
        self.query = query
        self.target = target

        # Alignment parameters
        self.ident = self.ins = self.dels = self.mis = 0
        self.score = score
        # Query and target lenghts.
        self.qlen = sum(1 for x in self.query.seq if x != '-')
        self.tlen = sum(1 for x in self.target.seq if x != '-')

        # Compute alignment information.
        for a, b in zip(self.query.seq, self.target.seq):
            if a == b:
                self.ident += 1
            elif a == '-':
                self.dels += 1
            elif b == '-':
                self.ins += 1
            elif a != b:
                self.mis += 1

        # Percent identity (BLAST identity)
        # https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
        self.pident = self.ident / (self.ident + self.mis + self.dels + self.ins) * 100

        self.alen = len(self.target.seq)


class Param:
    """
    Represents alignment parameters
    """

    def __init__(self, **kwds):
        self.is_dna = False
        self.matrix = None
        self.gap_open = 11
        self.gap_extend = 1
        self.mode = ''
        self.type = ''
        # Override with defaults
        self.__dict__.update(kwds)


MATCH, SNP, INS, DEL = 'M', 'SNP', 'INS', 'DEL'


def find_variants(query, target):
    """
    Takes two aligned sequences and tabulates what variants can be found at a given position.
    This is the trickiest to get right ... mostly works but may still have bugs.

    Returns a dictionary keyed by position where each position describes the  variant as a tuple.
    """

    stream = zip(count(), query.seq, target.seq, )

    # stream = islice(stream, 50000)

    collect = []
    variants = []
    lastop = None
    pos = 0
    for idx, a, b in stream:

        # Multiple sequence alignment
        if a == '-' and b == '-':
            continue

        if a == b:
            # Matches
            pos += 1
            op = MATCH

        elif a == '-':
            # Insertion into query
            pos += 1
            op = DEL

        elif b == '-':
            # Deletion from query.
            op = INS

        elif a != b:
            # Mismatching bases.
            pos += 1
            op = SNP

        else:
            raise Exception(f"Should never hit this (sanity check): {a} vs {b}")

        # Collect variants when the operator changes.
        if lastop != op:
            if collect:
                if lastop != MATCH:
                    variants.append((lastop, collect))
            # Reset the collector
            collect = []

        # Collect the positions that have been visited
        collect.append((idx, pos, a, b))

        # Keep track of the last previous operation
        lastop = op

    # Collect last element
    if collect:
        variants.append((lastop, collect))

    # This is necessary because consecutive variants may overlap SNP + INSERT for example
    vcfdict = dict()

    for key, elems in variants:

        if key == MATCH:
            continue

        # Lenght of variant
        size = len(elems)

        idx, pos = elems[0][0], elems[0][1]

        # Query starts with insertion
        pos = 1 if pos < 1 else pos

        ref = target.seq[idx:idx + size].strip('-')

        alt = query.seq[idx:idx + size].strip('-')

        ref = str(ref)
        alt = str(alt)

        info = f"TYPE={key}"

        if key == SNP:
            # Mismatches printed consecutively
            for idx, pos, alt, ref in elems:
                name = f"{pos}_{ref}_{alt}"
                value = [target.name, str(pos), name, ref, alt, ".", "PASS", info, "GT", "1"]
                vcfdict[pos] = value


        elif key == INS or key == DEL:
            # Handle insertions and deletions.

            # Push back on POS if it is not 1.
            if idx > 0:
                # For insertions the pos does not advance
                if key == DEL:
                    pos = pos - 1
                ref = target.seq[idx - 1] + ref
                alt = query.seq[idx - 1] + alt
            else:
                # When there is no preceding base the last base must be used
                lastidx = elems[-1][0] + 1
                ref = ref + target.seq[lastidx]
                alt = alt + query.seq[lastidx]

            alt = alt or '.'
            ref = ref or '.'

            if key == DEL:
                suff = len(ref[1:]) if len(ref) > 10 else ref[1:]
                name = f"{pos}_del_{suff}"
            else:
                suff = len(alt[1:]) if len(alt) > 10 else alt[1:]
                name = f"{pos}_ins_{suff}"

            value = [target.name, str(pos), name, ref, alt, ".", "PASS", info, "GT", "1"]
            vcfdict[pos] = value

    return vcfdict


def format_vcf(alns):
    for aln in alns:

        vcfdict = find_variants(query=aln.query, target=aln.target)

        query, target = aln.query, aln.target

        print('##fileformat=VCFv4.2')
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        print('##FILTER=<ID=PASS,Description="All filters passed">')
        print('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">')
        print(f'##contig=<ID={target.name},length={len(target.seq.strip("-"))},assembly={target.name}>')
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{query.name}")

        for value in vcfdict.values():
            print("\t".join(value))


def format_fasta(alns):
    for aln in alns:
        SeqIO.write([aln.query, aln.target], sys.stdout, "fasta")


def format_table(alns, sep="\t"):
    header = "target query pident len match mism ins del"
    print("\t".join(header.split()))

    for aln in alns:
        data = [
            f"{aln.target.name}", f"{aln.query.name}",
            f"{aln.pident:0.1f}", f"{aln.alen}",
            f"{aln.ident}", f"{aln.mis}", f"{aln.dels}", f"{aln.ins}"
        ]
        line = sep.join(data)
        print(line)


def format_diffs(alns):
    #header = "pos info query change target"
    #print("\t".join(header.split()))

    for idx, aln in enumerate(alns):

        vcfdict = find_variants(query=aln.query, target=aln.target)
        for value in vcfdict.values():
            name = value[2]
            pos = value[1]
            base = value[3]
            alt = value[4]
            size = '0'
            info = value[7]
            if SNP in info:
                info = SNP
                size = 1
            elif DEL in info:
                info = DEL
                base = base[1:] if pos != '1' else base[:1]
                alt = '-'
                size = len(base)

            elif INS in info:
                info = INS
                alt = alt[1:] if pos != '1' else alt[:1]
                base = '-'
                size = len(alt)

            data = [pos, info, aln.query.name, f"{alt}/{base}", aln.target.name, ]

            print("\t".join(data))


def format_pairwise(alns, par=None, width=81):
    """
    Formats an alignment in pairwise mode
    """
    for aln in alns:

        out = [
            "",
            f"# {aln.query.name} ({aln.qlen}) vs {aln.target.name} ({aln.tlen})",
            f"# pident={aln.pident:0.1f}% len={aln.alen} ident={aln.ident} mis={aln.mis} del={aln.dels} ins={aln.ins}",

        ]

        if par:
            score = f"score={aln.score} " if aln.score is not None else ''
            elem = f"# {par.mode}: {score}gap open={par.gap_open} extend={par.gap_extend}  matrix={par.matrix}"
            out.append(elem)

        out.append("")

        # Generate the traces
        for start in range(0, len(aln.target.seq), width):
            end = start + width

            seq1 = aln.query.seq[start:end]
            seq2 = aln.target.seq[start:end]

            trace = ""
            for a, b in zip(seq1, seq2):
                if a == b:
                    trace += '|'
                elif a == '-' or b == '-':
                    trace += '-'
                else:
                    trace += '.'

            out.append(str(seq1))
            out.append(trace)
            out.append(str(seq2))
            out.append("")

        print("\n".join(out))
