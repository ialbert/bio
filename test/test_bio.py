
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from generate import run


def test_1(capsys):
    cmd = "bio data --delete ncov,ratg13"
    run(cmd, capsys=capsys, fname=None)

def test_2(capsys):
    cmd = "bio fetch NC_045512 --rename ncov --seqid ncov"
    run(cmd, capsys=capsys, fname=None)

def test_3(capsys):
    cmd = "bio fetch MN996532 --rename ratg13 --seqid ratg13"
    run(cmd, capsys=capsys, fname=None)

def test_4(capsys):
    cmd = "bio convert ncov --json"
    run(cmd, capsys=capsys, fname="ncov.json")

def test_5(capsys):
    cmd = "bio convert ncov --genbank"
    run(cmd, capsys=capsys, fname="ncov.gb")

def test_6(capsys):
    cmd = "bio convert ncov --fasta"
    run(cmd, capsys=capsys, fname="ncov.fa")

def test_7(capsys):
    cmd = "bio convert ncov --gff"
    run(cmd, capsys=capsys, fname="ncov.gff")

def test_8(capsys):
    cmd = "bio convert ncov --gff --match phosphoesterase"
    run(cmd, capsys=capsys, fname="match.gff")

def test_9(capsys):
    cmd = "bio convert ncov --gff --gene S"
    run(cmd, capsys=capsys, fname="gene.gff")

def test_10(capsys):
    cmd = "bio convert ncov --gff  --start 10000 --end 20000"
    run(cmd, capsys=capsys, fname="overlap.gff")

def test_11(capsys):
    cmd = "bio convert ncov --gff  --start 10,000 --end 20,000"
    run(cmd, capsys=capsys, fname="overlap.gff")

def test_12(capsys):
    cmd = "bio convert ncov --gff  --start 10kb --end 20kb"
    run(cmd, capsys=capsys, fname="overlap.gff")

def test_13(capsys):
    cmd = "bio convert ncov --gff  --type CDS"
    run(cmd, capsys=capsys, fname="cds.gff")

def test_14(capsys):
    cmd = "bio convert ncov --gff  --type gene,CDS,mRNA"
    run(cmd, capsys=capsys, fname="manytypes.gff")

def test_15(capsys):
    cmd = "bio convert ncov --fasta --seqid foo --start 10 --end 20"
    run(cmd, capsys=capsys, fname="start.fa")

def test_16(capsys):
    cmd = "bio convert ncov --fasta --type CDS -end 10"
    run(cmd, capsys=capsys, fname="cds.fa")

def test_17(capsys):
    cmd = "bio convert ncov --fasta --type gene --end 10"
    run(cmd, capsys=capsys, fname="start-gene.fa")

def test_18(capsys):
    cmd = "bio convert ncov --fasta --translate --type CDS --end 10"
    run(cmd, capsys=capsys, fname="translates.fa")

def test_19(capsys):
    cmd = "bio convert ncov --fasta --protein --start -10"
    run(cmd, capsys=capsys, fname="protein-end.fa")

def test_20(capsys):
    cmd = "bio convert ncov --fasta --type CDS --gene S --end 10"
    run(cmd, capsys=capsys, fname="cds-gene.fa")

def test_21(capsys):
    cmd = "bio convert ncov:S --fasta --end 10"
    run(cmd, capsys=capsys, fname="cds-gene.fa")

def test_22(capsys):
    cmd = "bio convert ncov --id YP_009724390.1 --fasta --end 10"
    run(cmd, capsys=capsys, fname="cds-gene.fa")

def test_23(capsys):
    cmd = "bio convert ncov:S --fasta --protein --seqid foo"
    run(cmd, capsys=capsys, fname="cds-prot.fa")

def test_24(capsys):
    cmd = "bio convert ATGGGC -i --fasta"
    run(cmd, capsys=capsys, fname="inter.fa")

def test_25(capsys):
    cmd = "bio convert ATGGGC -i --translate --seqid foo"
    run(cmd, capsys=capsys, fname="inter-trans.fa")

def test_26(capsys):
    cmd = "bio convert ATGGGC -i --revcomp --translate --seqid foo"
    run(cmd, capsys=capsys, fname="inter-revcomp1.fa")

def test_27(capsys):
    cmd = "bio convert ATGGGC -i --reverse --complement --translate --seqid foo"
    run(cmd, capsys=capsys, fname="inter-revcomp2.fa")

def test_28(capsys):
    cmd = "bio align ncov ratg13 --end 180"
    run(cmd, capsys=capsys, fname="align-dna.txt")

def test_29(capsys):
    cmd = "bio align ncov:S ratg13:S --end 180"
    run(cmd, capsys=capsys, fname="align-gene.txt")

def test_30(capsys):
    cmd = "bio align ncov:S ratg13:S --end 180 --translate"
    run(cmd, capsys=capsys, fname="align-translation.txt")

def test_31(capsys):
    cmd = "bio align ratg13:S ncov:S  --end 180 -1"
    run(cmd, capsys=capsys, fname="align-gene-pept1.txt")

def test_32(capsys):
    cmd = "bio align ratg13:S ncov:S  --end 180 -3"
    run(cmd, capsys=capsys, fname="align-gene-pept3.txt")

def test_33(capsys):
    cmd = "bio align ncov:S ratg13:S --protein --end 60"
    run(cmd, capsys=capsys, fname="align-protein.txt")

def test_34(capsys):
    cmd = "bio align THISLINE ISALIGNED  -i --local"
    run(cmd, capsys=capsys, fname="align-local.txt")

def test_35(capsys):
    cmd = "bio align THISLINE ISALIGNED -i --global"
    run(cmd, capsys=capsys, fname="align-global.txt")

def test_36(capsys):
    cmd = "bio align THISLINE ISALIGNED -i --semiglobal"
    run(cmd, capsys=capsys, fname="align-semiglobal.txt")

def test_37(capsys):
    cmd = "bio taxon 9606"
    run(cmd, capsys=capsys, fname="taxon_9606.txt")

def test_38(capsys):
    cmd = "bio taxon 9606 --lineage"
    run(cmd, capsys=capsys, fname="taxon_9606_lineage.txt")

def test_39(capsys):
    cmd = "bio taxon 9606 --depth 1"
    run(cmd, capsys=capsys, fname="taxon_9606_depth_1.txt")

def test_40(capsys):
    cmd = "bio taxon  ncov ratg13"
    run(cmd, capsys=capsys, fname="taxon_data.txt")

def test_41(capsys):
    cmd = "bio data --delete ebola"
    run(cmd, capsys=capsys, fname=None)

def test_42(capsys):
    cmd = "bio fetch KM233118 --rename ebola"
    run(cmd, capsys=capsys, fname=None)

def test_43(capsys):
    cmd = "bio runinfo ebola"
    run(cmd, capsys=capsys, fname="sra-test.txt")

