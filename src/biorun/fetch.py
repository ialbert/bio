import csv
import os
import re
import sys
import tempfile
import zipfile
from pathlib import Path

try:
    from Bio import Entrez
except ImportError as exc:
    print(f"# Error: {exc}", file=sys.stderr)
    print(f"# This program requires biopython", file=sys.stderr)
    print(f"# Install: conda install -y biopython>=1.79", file=sys.stderr)
    sys.exit(-1)

from biorun import patterns
from biorun.api import ena_fastq

from biorun.libs import placlib as plac
from tqdm import tqdm
from biorun import utils
from urllib.error import HTTPError
import requests

Entrez.email = 'foo@foo.com'

logger = utils.logger


#
# Genbank and Refseq accession numbers
#
# https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
#

# https://rest.ensembl.org/sequence/id/ENST00000288602?content-type=text/x-fasta;type=cdna


def parse_ensmbl(text):
    patt = r'(?P<letters>[a-zA-Z]+)(?P<digits>\d+)(\.(?P<version>\d+))?'
    patt = re.compile(patt)
    m = patt.search(text)
    code = m.group("letters") if m else ''
    digits = m.group("digits") if m else ''
    version = m.group("version") if m else ''
    return code, digits, version


def parse_ncbi(text):
    patt = r'(?P<letters>[a-zA-Z]+)(?P<under>_?)(?P<digits>\d+)(\.(?P<version>\d+))?'
    patt = re.compile(patt)
    m = patt.search(text)
    code = m.group("letters") if m else ''
    digits = m.group("digits") if m else ''
    refseq = m.group("under") if m else ''
    version = m.group("version") if m else ''
    return code, digits, refseq, version


def is_GCF(text):
    return bool(re.fullmatch(r"GC[AF]_\d{9}(?:\.\d+)?", text, re.IGNORECASE))


def is_bioproject(text):
    """
    Pattern for project numbers.
    """
    return text.startswith("PRJ")


def is_ensembl(text):
    code, digits, version = parse_ensmbl(text)
    cond = code in ("ENST", "ENSG", "ENSP", "ENSE")
    cond = cond and len(digits) > 8 and digits.startswith('0')
    return cond


def is_ncbi_nucleotide(text):
    """
    Returns true of text matches NCBI nucleotides.
    """
    code, digits, refseq, version = parse_ncbi(text)
    if refseq:
        cond = code in ["AC", "NC", "NG", "NT", "NW", "NZ", "NM", "XM", "XR"]
    else:
        num1, num2 = len(code), len(digits)
        cond = (num1 == 1 and num2 == 5) or (num1 == 2 and num2 == 6) or (num1 == 3 and num2 == 8)
    return cond


def is_ncbi_protein(text):
    """
    Returns true of text matches NCBI protein sequences
    """
    code, digits, refseq, version = parse_ncbi(text)
    if refseq:
        cond = code in ["AP", "NP", "YP", "XP", "WP"]
    else:
        num1, num2 = len(code), len(digits)
        cond = (num1 == 3 and num2 == 5) or (num1 == 3 and num2 == 7)
    return cond


def efetch(ids, db, rettype, retmode, limit=None):
    """
    Fetches data from GenBank
    """
    ids = ",".join(ids) if type(ids) == list else ids

    try:
        stream = Entrez.efetch(db='sra', id=ids, rettype=rettype, retmode=retmode, retmax=limit)
    except HTTPError as exc:
        utils.error(f"HTTP Error: ids={ids}, db={db}, rettype={rettype}, retmode={retmode}, exc={exc}")
        stream = None

    return stream


def esearch_and_efetch(term, db, rettype, retmode, limit=None):
    """
    Search piped into a fetch
    """
    try:
        search = Entrez.esearch(db=db, term=term, usehistory="y", retmax=limit)
        results = Entrez.read(search)
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
        logger.debug(f"query={results}")
        stream = Entrez.efetch(db=db, webenv=webenv, query_key=query_key, retmax=limit, rettype=rettype,
                               retmode=retmode)
    except Exception as exc:
        utils.error(f"HTTP Error: term={term}, db={db}, rettype={rettype}, retmode={retmode}, exc={exc}")
        stream = None

    return stream


def format_results(stream, ftype=None):
    for line in stream:
        print(line, end='')
    return


def safe_int(text):
    try:
        return int(text)
    except ValueError:
        return 0


def format_runinfo(stream, ftype=None):
    if hasattr(stream, 'peek'):
        logger.debug(stream.read())
        logger.error("no valid results returned")
        sys.exit(1)

    # Get the CSV header.
    reader = csv.DictReader(stream)

    # Output may contain errors, empty rows and
    items = filter(lambda x: len(x.get("Run")) > 5, reader)

    if ftype == 'tsv':
        writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(items)
    if ftype == 'csv':
        writer = csv.DictWriter(sys.stdout, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(items)

    else:
        for elem in items:
            # print(elem)
            g = elem.get

            # Human friendly size
            size = safe_int(g('size_MB'))
            if size < 1024:
                size = f"{size:,d} MB"
            else:
                size = size / 1024
                size = f"{size:0.1f} GB"

            data = {
                "Project": f"{g('BioProject')}",
                "Run": g('Run'),
                "Library": f"{g('LibraryLayout')}, {g('LibrarySource')}, {g('LibraryStrategy')}",
                "Origin": f"{g('ScientificName')} ({g('TaxID')})",
                "Reads": f"{int(g('spots')):,d} (avgLength={g('avgLength')})",
                "Size": f"{size}",
                "Instr": f"{g('Platform')} ({g('Model')})",
                "Date": f"{g('LoadDate')}",
                "Path": f"{g('download_path')}",
            }
            print("")
            for key, value in data.items():
                print(f"{key}\t{value}")
        print("")


def fetch_ncbi_gff(ids, db="nuccore", chunk_size=1024*1024*10):
    url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi"
    params = dict(db=db, report='gff3', id=ids)
    try:
        r = requests.get(url, params=params)
        for chunk in r.iter_content(chunk_size=chunk_size):
            text = chunk.decode("utf8")
            print(text, sep='')
    except Exception as exc:
        utils.error(f"Error for {ids}, {db} : {exc}")


def fetch_ncbi(ids, db, rettype='gbwithparts', retmode='text', limit=None):
    ids = ",".join(ids) if type(ids) == list else ids

    try:
        stream = Entrez.efetch(db=db, id=ids, rettype=rettype, retmode=retmode, retmax=limit)
        bar = None
        use_bar = sys.stdout.isatty() and sys.stderr.isatty()

        if use_bar:
            try:
                bar = tqdm(stream, unit='B', unit_divisor=1024, desc='# downloaded', unit_scale=True, delay=5,
                           leave=False)
            except Exception as exc:
                # Older version of tqdm does not have the delay parameter.
                bar = tqdm(stream, unit='B', unit_divisor=1024, desc='# downloaded', unit_scale=True, leave=False)

    except Exception as exc:
        utils.error(f"Error for {ids}, {db}, {rettype}, {retmode}: {exc}")

    data = bar if bar else stream
    for line in data:
        print(line, end='')
        if bar:
            bar.update(len(line))

    if bar:
        bar.close()


DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"

ASSEMBLY_FILES = {
    "GENOME_FASTA": ("genomic.fna", "genome.fna"),
    "GENOME_GBFF": ("genomic.gbff", "genbank.gbff"),
    "GENOME_GFF": ("genomic.gff", "ann.gff"),
    "GENOME_GTF": ("genomic.gtf", "ann.gtf"),
    "CDS_FASTA": ("cds_from_genomic.fna", "cds.fna"),
    "RNA_FASTA": ("rna.fna", "rna.fna"),
    "PROT_FASTA": ("protein.faa", "protein.faa"),
    "SEQUENCE_REPORT": ("sequence_report.jsonl", "sequence_report.jsonl"),
}

DEFAULT_ASSEMBLY_INCLUDE = "genome,gff,gtf"

ASSEMBLY_INCLUDES = {
    "genome": "GENOME_FASTA",
    "gbff": "GENOME_GBFF",
    "genbank": "GENOME_GBFF",
    "gff": "GENOME_GFF",
    "gtf": "GENOME_GTF",
    "cds": "CDS_FASTA",
    "rna": "RNA_FASTA",
    "protein": "PROT_FASTA",
    "prot": "PROT_FASTA",
    "report": "SEQUENCE_REPORT",
    "sequence_report": "SEQUENCE_REPORT",
}

def assembly_types(include=DEFAULT_ASSEMBLY_INCLUDE):
    """Return NCBI Datasets annotation types selected by a comma-separated list."""
    names = [name.strip().lower() for name in include.split(",") if name.strip()]
    invalid = sorted(set(names) - set(ASSEMBLY_INCLUDES))
    if invalid:
        valid = ", ".join(sorted(ASSEMBLY_INCLUDES))
        utils.error(f"Invalid --include value: {', '.join(invalid)}. Valid values: {valid}")
    if not names:
        utils.error("--include requires at least one value")
    return {ASSEMBLY_INCLUDES[name] for name in names}


def format_mb(size):
    return f"{round(size / (1024 * 1024)):,}Mb"


def save_assembly(zip_path, acc, out, types, rename=""):
    """Copy selected assembly files from a Datasets ZIP to normalized paths."""
    names = {src: dest for kind, (src, dest) in ASSEMBLY_FILES.items() if kind in types}
    saved = []

    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            src = Path(member.filename).name
            suffix = names.get(src)
            if not suffix and re.fullmatch(
                    r"GC[AF]_\d{9}(?:\.\d+)?_.+_genomic\.fna", src, re.IGNORECASE):
                suffix = names.get("genomic.fna")
            if not suffix or member.is_dir():
                continue

            match = re.search(r"/(GC[AF]_\d{9}(?:\.\d+)?)/", member.filename, re.IGNORECASE)
            name = rename or (match.group(1).upper() if match else acc)
            dest = out / f"{name}.{suffix}"
            with archive.open(member) as inp, open(dest, "wb") as stream:
                while chunk := inp.read(1024 * 1024):
                    stream.write(chunk)
            saved.append(dest)

    return saved


def fetch_gcf(ids, include=DEFAULT_ASSEMBLY_INCLUDE, out="ncbi", rename=""):
    """Download GCF/GCA assemblies with the NCBI Datasets API."""
    if rename and len(ids) != 1:
        utils.error("--rename requires exactly one GCF/GCA accession")
    if rename and Path(rename).name != rename:
        utils.error("--rename must be a file name, not a path")

    out = Path(out or "ncbi")
    out.mkdir(parents=True, exist_ok=True)
    types = assembly_types(include)
    params = [("include_annotation_type", value) for value in sorted(types)]
    downloads = []

    for value in ids:
        acc = value.upper()
        url = f"{DATASETS_URL}/{acc}/download"
        temp = None
        try:
            with requests.get(url, params=params, headers={"Accept": "application/zip"},
                              stream=True, timeout=120) as response:
                response.raise_for_status()
                size = int(response.headers.get("content-length", 0))
                bar = tqdm(total=size or None, unit="B", unit_scale=True, unit_divisor=1024,
                           desc=f"# downloading {acc}")
                with tempfile.NamedTemporaryFile(delete=False, suffix=".zip") as stream, bar:
                    temp = stream.name
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            stream.write(chunk)
                            bar.update(len(chunk))

            saved = save_assembly(temp, acc, out, types, rename=rename)
            if not saved:
                utils.error(f"No requested files found for {acc}")
            downloads.extend(saved)
        except (requests.RequestException, OSError, zipfile.BadZipFile) as exc:
            utils.error(f"Error downloading {acc}: {exc}")
        finally:
            if temp and os.path.exists(temp):
                os.unlink(temp)

    for path in downloads:
        print(f"# {path} ({format_mb(path.stat().st_size)})")


def fetch_ensembl(ids, ftype='genomic'):
    ftype = 'genomic' if not ftype else ftype

    server = "https://rest.ensembl.org"

    try:
        for acc in ids:

            ext = f"/sequence/id/{acc}?type={ftype}"

            r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            print(r.text)
    except Exception as exc:
        utils.error(f"Error for {ids}, {ftype}: {exc}")


@plac.pos("acc", "accession numbers")
@plac.opt("db", "database", choices=["nuccore", "protein"])
@plac.opt("format_", "return format", choices=["gbwithparts", "fasta", "gb", "csv", "tsv", "default",
                                               "gff", "gtf", "cds", "rna", "prot", "protein"])
@plac.opt("type_", "get CDS/CDNA (Ensembl only)")
@plac.opt("limit", "limit results")
@plac.opt("out", "output directory for GCF/GCA files (default: ncbi)")
@plac.opt("rename", "rename the GCF/GCA output prefix")
@plac.opt("include", "GCF/GCA files: genome,gff,gtf,gbff,cds,rna,protein,report")
def run(db="", type_='', format_='', limit=100, out='', rename='',
        include=DEFAULT_ASSEMBLY_INCLUDE, *acc):
    ids = []
    for num in acc:
        ids.extend(num.split(","))

    if not sys.stdin.isatty():
        lines = utils.read_lines(sys.stdin, sep=None)
        ids.extend(lines)

    srrs = list(map(patterns.is_srr, ids))
    if all(srrs):
        for srr in ids:
            ena_fastq.run(srr=srr, out=out, limit=limit)
        return

    # Dealing with Ensembl
    ensmbl = list(map(is_ensembl, ids))
    if all(ensmbl):
        fetch_ensembl(ids=ids, ftype=type_)
        return

    # Downloads GCF data
    gcfs = list(map(is_GCF, ids))
    if all(gcfs):
        fetch_gcf(ids=ids, include=include, out=out or "ncbi", rename=rename)
        return

    # Detects nucleotides
    nucs = list(map(is_ncbi_nucleotide, ids))

    # Detects proteins
    prots = list(map(is_ncbi_protein, ids))

    if any(prots) and any(nucs):
        utils.error(f"input mixes protein and nucleotide entries: {ids}")

    # Select protein database
    if all(prots):
        default = "protein"
    else:
        default = "nuccore"

    db = default if not db else db
    # Fetch the ids
    ids = ",".join(ids)

    # Set a default.
    format_ = "gbwithparts" if not format_ else format_

    if format_ == "gff":
        fetch_ncbi_gff(ids=ids)
    elif format_ in ("genbank", "fasta", "gbwithparts"):
        rettype = "gbwithparts" if format_ == "genbank" else format_
        fetch_ncbi(db=db, rettype=rettype, ids=ids)
    else:
        utils.error(f"Not a valid format for this accession number: {format_}")


if __name__ == '__main__':
    # id = "AY851612",

    run()
