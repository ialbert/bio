import sys
from itertools import *

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def positions(query, target):
    pos = dict()

    posa = posb = 0

    for a, b in zip(query, target):

        if a != b:
            pos[posb] = posa

        if a != '-':
            posa += 1

        if b != '-':
            posb += 1

    return pos


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

        # Query and target lengths.
        self.qlen = sum(1 for x in self.query.seq if x != '-')
        self.tlen = sum(1 for x in self.target.seq if x != '-')

        # Compute alignment information.
        for a, b in zip(self.query, self.target):
            if a == '-' and b == '-':
                continue
            elif a == b:
                self.ident += 1
            elif a == '-':
                self.dels += 1
            elif b == '-':
                self.ins += 1
            elif a != b:
                self.mis += 1

        # Maps the positions between the two sequences
        self.pos = positions(query=self.query, target=self.target)

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
        self.match = 1
        self.mismatch = 2
        self.gap_open = 11
        self.gap_extend = 1
        self.is_dna = True
        self.mode = ''
        self.type = ''
        # Override with defaults
        self.__dict__.update(kwds)


MATCH, SNP, INS, DEL, VAR, SUB = 'MTC', 'SNP', 'INS', 'DEL', 'VAR', 'SUB'


def find_variants(query, target):
    """
    Takes two aligned sequences and tabulates what variants can be found at a given position.
    This is the trickiest to get right ... mostly works but may still have bugs.

    Returns a dictionary keyed by position where each position describes the  variant as a tuple.
    """

    stream = zip(count(), query.seq, target.seq, )

    # stream = islice(stream, 50000)

    vara = varb = ''
    variants = []
    lastop = None
    pos = 0
    for idx, a, b in stream:

        # Multiple sequence alignment
        if a == '-' and b == '-':
            continue

        if a == b:
            # Matches
            op = MATCH

        elif a == '-':
            # Deletion from query.
            op = DEL

        elif b == '-':
            # Insertion into query.
            op = INS

        elif a != b:
            # Mismatching bases.
            op = SNP
        else:
            raise Exception(f"Should never hit this (sanity check): {a} vs {b}")

        # Collect variants when the operator changes.
        if lastop != op:
            # Collect subsequences when we see a match only
            if op == MATCH or lastop == MATCH:
                if vara:
                    label = lastop if lastop == MATCH else VAR
                    variants.append((label, pos, vara, varb))
                    pos += len(varb) - varb.count('-')
                vara = varb = ''

        # Keep track of the last previous operation
        lastop = op
        vara += a
        varb += b

    # Collect last element
    if vara:
        variants.append((lastop, pos, vara, varb))

    last_a, last_b = None, None

    store = []
    for key, pos, seqa, seqb in variants:

        if key == MATCH:
            info = f"TYPE={MATCH}"
            name = f"{pos + 1}_{seqb}/{seqa}"
            row = [target.name, str(pos + 1), name, seqb, seqa, "PASS", info, "GT", "1"]
            last_a = seqa[-1]
            last_b = seqb[-1]
            # store.append(row)
            continue

        gaps_a = seqa.count("-")
        gaps_b = seqb.count("-")

        if not gaps_a and not gaps_b:
            # SNPS
            for i, a, b in zip(count(), seqa, seqb):
                info = f"TYPE={SNP}"
                name = f"{pos + i + 1}_{SNP}_{b}/{a}"
                row = [target.name, str(pos + i + 1), name, b, a, ".", "PASS", info, "GT", "1"]
                store.append(row)
            continue

        alt = seqa.replace("-", "")
        ref = seqb.replace("-", "")

        if alt == ref:
            # Rarely but happens.
            # Alignment gone bad.
            continue

        if last_b:
            alt = last_a + alt
            ref = last_b + ref
        else:
            pos = 1
            str1 = str(query.seq).replace("-", "")
            str2 = str(target.seq).replace("-", "")
            alt = alt + str1[len(alt)]
            ref = ref + str2[len(ref)]

        if len(ref) == 1:
            # Insertion
            info = f"TYPE={INS}"
            name = f"{pos}_{INS}_{len(alt) - 1}"

        elif len(alt) == 1:
            # Deletion
            info = f"TYPE={DEL}"
            name = f"{pos}_{DEL}_{len(ref) - 1}"

        else:
            info = f"TYPE={SUB}"
            name = f"{pos}_{SUB}_{len(ref) - 1}"

        row = [target.name, str(pos), name, ref, alt, ".", "PASS", info, "GT", "1"]

        store.append(row)

    return store


def format_vcf(alns):
    for aln in alns:

        values = find_variants(query=aln.query, target=aln.target)

        query, target = aln.query, aln.target

        print('##fileformat=VCFv4.2')
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        print('##FILTER=<ID=PASS,Description="All filters passed">')
        print('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of the variant">')
        print(f'##contig=<ID={target.name},length={len(target.seq.strip("-"))},assembly={target.name}>')
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{query.name}")

        for row in values:
            print("\t".join(row))


def format_fasta(alns):
    for aln in alns:
        SeqIO.write([aln.target, aln.query], sys.stdout, "fasta")


def format_table(alns, sep="\t"):
    header = "target query pident len match mism ins del"

    print("\t".join(header.split()))

    for aln in alns:
        data = [
            f"{aln.target.name}", f"{aln.query.name}",
            f"{aln.pident:0.1f}", f"{aln.alen}",
            f"{aln.ident}", f"{aln.mis}", f"{aln.ins}", f"{aln.dels}"
        ]
        line = sep.join(data)
        print(line)


def format_mutations(alns):
    # header = "pos info query change target"
    # print("\t".join(header.split()))

    for idx, aln in enumerate(alns):

        values = find_variants(query=aln.query, target=aln.target)

        #print(f"# {aln.target.name} vs {aln.query.name}")

        for elems in values:
            pos = elems[1]
            base = elems[3]
            alt = elems[4]
            info = elems[7].split("=")[-1]
            short = f"{base}{pos}{alt}"

            # data = [short, pos, info, aln.target.name, f"{base}/{alt}", aln.query.name, ]

            data = [short, info, pos, base, alt, aln.target.name, aln.query.name]

            print("\t".join(data))


def format_pile(alns):
    coords = dict()
    data = dict()
    names = dict()

    for idx, aln in enumerate(alns):

        names[aln.target.name] = True
        names[aln.query.name] = True

        values = find_variants(query=aln.query, target=aln.target)

        data[idx] = dict()

        for elems in values:

            pos = elems[1]
            base = elems[3]
            alt = elems[4]
            info = elems[7].split("=")[-1]
            short = f"{base}{pos}{alt}"
            pos = int(pos)
            if info != SNP:
                continue

            coords[(pos, base, short)] = True

            data[idx][short] = alt

    header = " vs ".join(names.keys())
    print(f"# {header}")
    for coord in sorted(coords.keys()):

        pos, base, short = coord

        out = [short, f"{pos}", f"{base}"]

        found = 0
        pile = ''
        for idx in data:
            alt = data[idx].get(short, ".")
            if alt != '.':
                found += 1

            pile += alt

        out.append(pile)
        out.append(f"{found}")

        print("\t".join(out))


def format_pairwise(alns, par=None, width=81):
    """
    Formats an alignment in pairwise mode
    """
    for aln in alns:

        out = [
            "",
            f"# {aln.target.name} ({aln.tlen}) vs {aln.query.name} ({aln.qlen})",
            f"# pident={aln.pident:0.1f}% len={aln.alen} ident={aln.ident} mis={aln.mis} del={aln.dels} ins={aln.ins}",
        ]

        if par:
            score = f"score={aln.score}" if aln.score is not None else ''
            if par.matrix:
                elem = f"# {par.mode}: {score} matrix={par.matrix} gapopen={par.gap_open} gapextend={par.gap_extend}"
            else:
                elem = f"# {par.mode}: {score} match={par.match} mismatch={par.mismatch} gapopen={par.gap_open} gapextend={par.gap_extend}"
            out.append(elem)

        out.append("")

        # Generate the traces
        for start in range(0, len(aln.target.seq), width):
            end = start + width

            seq1 = aln.target.seq[start:end]
            seq2 = aln.query.seq[start:end]

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
