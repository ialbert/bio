
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from generate import run


def test_1(capsys):
    cmd = "bio ncov --delete"
    run(cmd, capsys=capsys, fname=None)

def test_2(capsys):
    cmd = "bio NC_045512 --fetch --rename ncov --seqid ncov"
    run(cmd, capsys=capsys, fname=None)

def test_3(capsys):
    cmd = "bio NC_045512 --fasta --end 100"
    run(cmd, capsys=capsys, fname="fetch-convert.fa")

def test_4(capsys):
    cmd = "bio ncov"
    run(cmd, capsys=capsys, fname="ncov.json")

def test_5(capsys):
    cmd = "bio ncov --genbank"
    run(cmd, capsys=capsys, fname="ncov.gb")

def test_6(capsys):
    cmd = "bio ncov --genome"
    run(cmd, capsys=capsys, fname="ncov-genome.fa")

def test_7(capsys):
    cmd = "bio ncov --fasta"
    run(cmd, capsys=capsys, fname="ncov-features.fa")

def test_8(capsys):
    cmd = "bio ncov --gff"
    run(cmd, capsys=capsys, fname="ncov.gff")

def test_9(capsys):
    cmd = "bio ncov --gff --match ORF1ab"
    run(cmd, capsys=capsys, fname="match.gff")

def test_10(capsys):
    cmd = "bio ncov --gff --gene S"
    run(cmd, capsys=capsys, fname="gene1.gff")

def test_11(capsys):
    cmd = "bio ncov --gff  --start 10000 --end 20000"
    run(cmd, capsys=capsys, fname="overlap.gff")

def test_12(capsys):
    cmd = "bio ncov --gff  --start 10,000 --end 20,000"
    run(cmd, capsys=capsys, fname="overlap.gff")

def test_13(capsys):
    cmd = "bio ncov --gff  --start 10kb --end 20kb"
    run(cmd, capsys=capsys, fname="overlap.gff")

def test_14(capsys):
    cmd = "bio ncov --gff  --type CDS"
    run(cmd, capsys=capsys, fname="type.gff")

def test_15(capsys):
    cmd = "bio ncov --gff  --type gene,CDS,mRNA"
    run(cmd, capsys=capsys, fname="multiple-types.gff")

def test_16(capsys):
    cmd = "bio ncov --fasta --seqid foo --start 10 --end 20"
    run(cmd, capsys=capsys, fname="fasta-start.fa")

def test_17(capsys):
    cmd = "bio ncov --fasta --type CDS"
    run(cmd, capsys=capsys, fname="CDS.fa")

def test_18(capsys):
    cmd = "bio ncov --fasta --type gene --end 10"
    run(cmd, capsys=capsys, fname="gene-start.fa")

def test_19(capsys):
    cmd = "bio ncov --fasta --translate --type CDS"
    run(cmd, capsys=capsys, fname="translate.fa")

def test_20(capsys):
    cmd = "bio ncov --fasta --protein --start -10"
    run(cmd, capsys=capsys, fname="protein-end.fa")

def test_21(capsys):
    cmd = "bio ncov --fasta --type CDS --gene S --end 10"
    run(cmd, capsys=capsys, fname="cds-gene-s.fa")

def test_22(capsys):
    cmd = "bio ncov:S --fasta --end 10"
    run(cmd, capsys=capsys, fname="cds-gene-s.fa")

def test_23(capsys):
    cmd = "bio ncov --id YP_009724390.1 --fasta --end 10"
    run(cmd, capsys=capsys, fname="cds-gene-s.fa")

def test_24(capsys):
    cmd = "bio ncov:S --fasta --protein --seqid foo"
    run(cmd, capsys=capsys, fname="s_prot_foo.fa")

def test_25(capsys):
    cmd = "bio ATGGGC -i --fasta"
    run(cmd, capsys=capsys, fname="inter.fa")

def test_26(capsys):
    cmd = "bio ATGGGC -i --translate --seqid foo"
    run(cmd, capsys=capsys, fname="inter-trans.fa")

def test_27(capsys):
    cmd = "bio ATGGGC -i --revcomp --translate --seqid foo"
    run(cmd, capsys=capsys, fname="inter-revcomp1.fa")

def test_28(capsys):
    cmd = "bio ATGGGC -i --reverse --complement --translate --seqid foo"
    run(cmd, capsys=capsys, fname="inter-revcomp2.fa")

def test_29(capsys):
    cmd = "bio MN996532 --fetch --rename ratg13 --seqid ratg13"
    run(cmd, capsys=capsys, fname=None)

def test_30(capsys):
    cmd = "bio ncov ratg13 --end 210 --align"
    run(cmd, capsys=capsys, fname="align-dna.txt")

def test_31(capsys):
    cmd = "bio ncov:S ratg13:S --end 210 --align"
    run(cmd, capsys=capsys, fname="align-dna-s.txt")

def test_32(capsys):
    cmd = "bio ncov:S ratg13:S --end 210 --translate --align"
    run(cmd, capsys=capsys, fname="align-translated-s.txt")

def test_33(capsys):
    cmd = "bio ncov:S ratg13:S --protein --end 70 --align"
    run(cmd, capsys=capsys, fname="align-protein-s.txt")

def test_34(capsys):
    cmd = "bio ratg13:S ncov:S  --start 91 --end 120 --align -1"
    run(cmd, capsys=capsys, fname="align-short-pept.txt")

def test_35(capsys):
    cmd = "bio ratg13:S ncov:S  --start 91 --end 120 --align -3"
    run(cmd, capsys=capsys, fname="align-long-pept.txt")

def test_36(capsys):
    cmd = "bio THISLINE ISALIGNED  -i --align --local"
    run(cmd, capsys=capsys, fname="align-local.txt")

def test_37(capsys):
    cmd = "bio THISLINE ISALIGNED -i --align --global"
    run(cmd, capsys=capsys, fname="align-global.txt")

def test_38(capsys):
    cmd = "bio THISLINE ISALIGNED -i --align --semiglobal"
    run(cmd, capsys=capsys, fname="align-semiglobal.txt")

def test_39(capsys):
    cmd = "bio 9606 --taxon"
    run(cmd, capsys=capsys, fname="taxon_9606.txt")

def test_40(capsys):
    cmd = "bio 9606 --lineage --taxon"
    run(cmd, capsys=capsys, fname="taxon_9606_lineage.txt")

def test_41(capsys):
    cmd = "bio 9606 --lineage --flat --taxon"
    run(cmd, capsys=capsys, fname="taxon_9606_flat_lineage.txt")

def test_42(capsys):
    cmd = "bio ncov ratg13 --taxon"
    run(cmd, capsys=capsys, fname="taxon_ncov_ratg13.txt")

def test_43(capsys):
    cmd = "bio ebola --delete"
    run(cmd, capsys=capsys, fname=None)

def test_44(capsys):
    cmd = "bio KM233118 --fetch --rename ebola"
    run(cmd, capsys=capsys, fname=None)

def test_45(capsys):
    cmd = "bio ebola --sra"
    run(cmd, capsys=capsys, fname="sra-test.txt")

