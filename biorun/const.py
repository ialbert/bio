# Constants reused in the program

# JSON field names
DEFINITION = "definition"
DBLINK = "dblink"
LOCUS = "locus"
SEQID = "id"
ORIGIN = "ORIGIN"
FEATURES = "FEATURES"
FEATURE_COUNT = "feature_count"
ORIGIN_SIZE = "origin_len"

# Alignment modes.
GLOBAL_ALIGN, LOCAL_ALIGN, SEMIGLOBAL_ALIGN, STRICT_ALIGN = "global", "local", "semiglobal", "strictglobal"

ALIGN, TAXON = "align", "taxon"

# Actions that the bio package recognizes
ACTIONS_WORDS = [ ALIGN ]

#
# Remaps types from GenBank to Sequence Ontology when converting to GFF files
#
SEQUENCE_ONTOLOGY = {
    "5'UTR": "five_prime_UTR",
    "mat_peptide": "mature_protein_region",
}

#
# The GFF attributes generated for a source type.
#
SOURCE_ATTRIBUTES = [
    "mol_type", "isolate", "db_xref", "organism", "country", "collection_date"
]

# GFF attributes filled for each feature other than "source"
GFF_ATTRIBUTES = [
    "gene", "protein_id", "product", "db_xref", "function",
]

# Recognized types.
GENBANK, FASTA, GFF, BED, SAM, BAM = "genbank", "fasta", "gff", "bed", "sam", "bam"

# Connect an extension to types.
TYPE_BY_EXTENSION = {
    "gb": GENBANK,
    "gbk": GENBANK,
    "genbank": GENBANK,
    "fa": FASTA,
    "fasta": FASTA,
    "bed": BED,
    "gff": GFF,
    "sam": SAM,
    "bam": BAM,
}