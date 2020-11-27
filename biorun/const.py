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
GLOBAL_ALIGN, LOCAL_ALIGN, SEMIGLOBAL_ALIGN, STRICT_GLOBAL_ALIGN = "global", "local", "semiglobal", "strictglobal"

ALIGN_COMMAND, TAXON_COMMAND, DBLINK_COMMAND = "--align", "--taxon", "--sra"

SKIP_GFF_ATTR = { "id", "parent_id", "name",  "type", "start", "end", "location", "translation", "strand", "operator"}


# Guess accession numbers that are proteins based on start letters
# https: // www.ncbi.nlm.nih.gov / Sequin / acc.html

#
# Remaps types from GenBank to Sequence Ontology when converting to GFF files
#

NCBI_PROTEIN_CODES = {"AP", "NP", "YP", "XP", "WP", "AK"}

NCBI_NUCLEOTIDE_CODES = {"NM", "NX", "YP", "XP", "WP", "AK"}


SEQUENCE_ONTOLOGY = {
    "source": "region",
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}

# Feature types where the name of the feature is predetermined.
NAME_FROM_TYPE = {
    "five_prime_UTR": "five_prime_UTR",
    "three_prime_UTR": "three_prime_UTR",
    "stem_loop": "stem_loop",
}

# Associates a color to a feature type.
COLOR_FOR_TYPE = {
    "five_prime_UTR": "#cc0e74",
    "three_prime_UTR": "#cc0e74",
    "stem_loop": "#fa7f72",
    "mature_protein_region": "#CBAEBB",
    "region": "#CECECE",
    "transcript": "#799351",
    "gene": "#cb7a77",
    "mRNA": "#7a77cb",
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