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

# Indentation character
INDENT = '  '

# Fields separator

SEP = ', '

CHUNK = 25000

# Alignment modes.
GLOBAL_ALIGN, LOCAL_ALIGN, SEMIGLOBAL_ALIGN, STRICT_GLOBAL_ALIGN = "global", "local", "semiglobal", "strictglobal"

# Command map
# module.function, helpflag, descriptin
SUB_COMMANDS = dict(
    data=("biorun.fetch.data", False, "list or rename data"),
    fetch=("biorun.fetch.run", True, "downloads data from repositories"),
    align=("biorun.methods.align.run", True, "performs sequence alignments"),
    taxon=("biorun.models.taxdb.run", False, "displays NCBI taxonomies"),
    define=("biorun.models.ontology.run", False, "explains biological terms"),
    convert=("biorun.convert.run", True, "converts data to different formats"),
    runinfo=("biorun.runinfo.run", True, "prints sequencing run information"),
)


SKIP_GFF_ATTR = {"id", "parent_id", "name", "type", "start", "end", "location", "translation", "strand", "operator"}

# Guess accession numbers that are proteins based on start letters
# https: // www.ncbi.nlm.nih.gov / Sequin / acc.html

#
# Remaps types from GenBank to Sequence Ontology when converting to GFF files
#

NCBI_PROTEIN_CODES = {"AP", "NP", "YP", "XP", "WP", "AK"}

NCBI_NUCLEOTIDE_CODES = {"NM", "NX", "YP", "XP", "WP", "AK"}

BUCKET_NAME = "biostore-bucket-001"

# Types with hierachies
MULTIPART_TYPES = {"mRNA", "CDS", "ncRNA", "tRNA"}

SEQUENCE_ONTOLOGY = {
    "source": "region",
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "mat_peptide": "mature_protein_region",
}

# Column headers for the genome assembly data.
GENOME_ASSEMBLY_HEADER = ['assembly_accession', 'bioproject', 'biosample',
                  'wgs_master', 'refseq_category', 'taxid',
                  'species_taxid', 'organism_name', 'infraspecific_name',
                  'isolate', 'version_status', 'assembly_level', 'release_type',
                  'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
                  'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
                  'excluded_from_refseq', 'relation_to_type_material'
                          ]

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
    "mRNA": "#799351",
    "gene": "#cb7a77",
    "transcript": "#79a3b1",
    "tRNA": "#a685e2",
    "ncRNA": "#fca3cc",
    "mobile_element": "#efd9d1",
    "mRNA_region":"#7a77cb",
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

# Map and association file to a given field
ASSOCIATION_MAP = {'aspgd': 'http://current.geneontology.org/annotations/aspgd.gaf.gz',
                  'cgd': 'http://current.geneontology.org/annotations/cgd.gaf.gz',
                  'dictybase': 'http://current.geneontology.org/annotations/dictybase.gaf.gz',
                  'ecocyc': 'http://current.geneontology.org/annotations/ecocyc.gaf.gz',
                  'fb': 'http://current.geneontology.org/annotations/fb.gaf.gz',
                  'genedb_lmajor': 'http://current.geneontology.org/annotations/genedb_lmajor.gaf.gz',
                   'genedb_tbrucei': 'http://current.geneontology.org/annotations/genedb_tbrucei.gaf.gz',
                   'chicken': 'http://current.geneontology.org/annotations/goa_chicken.gaf.gz',
                   'chicken_complex': 'http://current.geneontology.org/annotations/goa_chicken_complex.gaf.gz',
                   'chicken_isoform': 'http://current.geneontology.org/annotations/goa_chicken_isoform.gaf.gz',
                   'chicken_rna': 'http://current.geneontology.org/annotations/goa_chicken_rna.gaf.gz',
                   'cow': 'http://current.geneontology.org/annotations/goa_cow.gaf.gz',
                   'cow_complex': 'http://current.geneontology.org/annotations/goa_cow_complex.gaf.gz',
                   'cow_isoform': 'http://current.geneontology.org/annotations/goa_cow_isoform.gaf.gz',
                   'cow_rna': 'http://current.geneontology.org/annotations/goa_cow_rna.gaf.gz',
                   'dog': 'http://current.geneontology.org/annotations/goa_dog.gaf.gz',
                   'dog_complex': 'http://current.geneontology.org/annotations/goa_dog_complex.gaf.gz',
                   'dog_isoform': 'http://current.geneontology.org/annotations/goa_dog_isoform.gaf.gz',
                   'dog_rna': 'http://current.geneontology.org/annotations/goa_dog_rna.gaf.gz',
                   'human': 'http://current.geneontology.org/annotations/goa_human.gaf.gz',
                   'human_complex': 'http://current.geneontology.org/annotations/goa_human_complex.gaf.gz',
                   'human_isoform': 'http://current.geneontology.org/annotations/goa_human_isoform.gaf.gz',
                   'human_rna': 'http://current.geneontology.org/annotations/goa_human_rna.gaf.gz',
                   'pig': 'http://current.geneontology.org/annotations/goa_pig.gaf.gz',
                   'pig_complex': 'http://current.geneontology.org/annotations/goa_pig_complex.gaf.gz',
                   'pig_isoform': 'http://current.geneontology.org/annotations/goa_pig_isoform.gaf.gz',
                   'pig_rna': 'http://current.geneontology.org/annotations/goa_pig_rna.gaf.gz',
                   'uniprot_all': 'http://current.geneontology.org/annotations/goa_uniprot_all.gaf.gz',
                   'uniprot_all_noiea': 'http://current.geneontology.org/annotations/goa_uniprot_all_noiea.gaf.gz',
                   'mgi': 'http://current.geneontology.org/annotations/mgi.gaf.gz',
                   'pombase': 'http://current.geneontology.org/annotations/pombase.gaf.gz',
                   'pseudocap': 'http://current.geneontology.org/annotations/pseudocap.gaf.gz',
                   'reactome': 'http://current.geneontology.org/annotations/reactome.gaf.gz',
                   'rgd': 'http://current.geneontology.org/annotations/rgd.gaf.gz',
                   'sgd': 'http://current.geneontology.org/annotations/sgd.gaf.gz',
                   'sgn': 'http://current.geneontology.org/annotations/sgn.gaf.gz',
                   'tair': 'http://current.geneontology.org/annotations/tair.gaf.gz',
                   'wb': 'http://current.geneontology.org/annotations/wb.gaf.gz',
                   'zfin': 'http://current.geneontology.org/annotations/zfin.gaf.gz'
                   }