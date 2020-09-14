from django.db import models


MAX_CHAR = 256
MAX_TEXT = 100000


class Locus(models.Model):
    # https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

    # First line in file
    name = models.CharField(max_length=MAX_CHAR)
    definition = models.TextField(max_length=MAX_TEXT)

    # Accession number of this specific locus.
    accession = models.CharField(max_length=MAX_CHAR)

    version = models.CharField(max_length=MAX_CHAR)
    dblink = models.CharField(max_length=MAX_CHAR)
    keywords = models.CharField(max_length=MAX_CHAR)

    organism = models.TextField(max_length=MAX_TEXT)

    # Store orin as fasta file path.
    origin = models.FilePathField()

    comment = models.TextField(max_length=MAX_TEXT)

    def __str__(self):
        return self.name


class Feature(models.Model):

    # source, CDS, mRNA, gene, etc...
    type = models.CharField(max_length=MAX_CHAR)

    # Range of integers
    location = models.CharField(max_length=MAX_CHAR)
    start = ''
    end = ''

    gene = models.CharField(max_length=MAX_CHAR)

    locus_tag = models.CharField(max_length=MAX_CHAR)
    ribosomal_slippage = models.BooleanField(default=False)

    note = models.TextField(max_length=MAX_TEXT)

    codon_start = models.BooleanField(default=False)
    product = models.CharField(max_length=MAX_CHAR)
    protein_id = models.CharField(max_length=MAX_CHAR)

    db_xref = models.CharField(max_length=MAX_CHAR)
    inference = models.CharField(max_length=MAX_CHAR)
    function = models.CharField(max_length=MAX_CHAR)

    bound_moiety = models.CharField(max_length=MAX_CHAR)
    experiment = models.CharField(max_length=MAX_CHAR)
    translation = models.CharField(max_length=MAX_CHAR)

    EC_number = models.CharField(max_length=MAX_CHAR)
    direction = models.CharField(max_length=MAX_CHAR)

    operon = models.CharField(max_length=MAX_CHAR)
    number = models.IntegerField()
    clone = models.BooleanField(default=False)
    anticodon = models.CharField(max_length=MAX_CHAR)

    # Attributes specific to type 'repeat_region'
    rpt_family = models.CharField(max_length=MAX_CHAR)
    rpt_type = models.CharField(max_length=MAX_CHAR)
    rpt_unit_seq = models.CharField(max_length=MAX_CHAR)
    regulatory_class = models.CharField(max_length=MAX_CHAR)

    # Attributes specific to type 'source'
    organism = models.CharField(max_length=MAX_CHAR)
    mol_type = models.CharField(max_length=MAX_CHAR)
    ecotype = models.CharField(max_length=MAX_CHAR)
    strain = models.CharField(max_length=MAX_CHAR)
    chromosome = models.IntegerField()
    map = models.CharField(max_length=MAX_CHAR)
    lab_host = models.CharField(max_length=MAX_CHAR)
    dev_stage = models.CharField(max_length=MAX_CHAR)
    focus = models.BooleanField(default=False)

    rearranged = models.BooleanField(default=False)
    cell_line = models.CharField(max_length=MAX_CHAR)
    cell_type = models.CharField(max_length=MAX_CHAR)

    # Get features associated with locus
    locus = models.ForeignKey(Locus, on_delete=models.SET_NULL, null=True)

    def sources(self):
        return

    def CDS(self):
        return

    def regulatory(self):
        return

    def repeat_region(self):
        return







