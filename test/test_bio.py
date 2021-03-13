
#
# This file was generated automatically! Do not edit.
#

# Get the helper utitilies.
from generate import *

# Initialize directories.
init_dirs()


def test_1(capsys):
    run("set -uex")

def test_2(capsys):
    run("bio data --delete ncov,ratg13")

def test_3(capsys):
    run("bio fetch NC_045512 --rename ncov --seqid ncov")

def test_4(capsys):
    run("bio fetch MN996532 --rename ratg13 --seqid ratg13")

def test_5(capsys):
    run("bio convert ncov --json > ncov.json")

def test_6(capsys):
    run("bio convert ncov --genbank > ncov.gb")

def test_7(capsys):
    run("bio convert ncov --fasta > ncov.fa")

def test_8(capsys):
    run("bio convert ncov --gff > ncov.gff")

def test_9(capsys):
    run("bio convert ncov --gff --match phosphoesterase  > match.gff")

def test_10(capsys):
    run("bio convert ncov --gff --gene S > gene.gff")

def test_11(capsys):
    run("bio convert ncov --gff  --start 10000 --end 20000 > overlap.gff")

def test_12(capsys):
    run("bio convert ncov --gff  --start 10,000 --end 20,000 > overlap.gff")

def test_13(capsys):
    run("bio convert ncov --gff  --start 10kb --end 20kb > overlap.gff")

def test_14(capsys):
    run("bio convert ncov --gff  --type CDS > cds.gff")

def test_15(capsys):
    run("bio convert ncov --gff  --type gene,CDS,mRNA > manytypes.gff")

def test_16(capsys):
    run("bio convert ncov --fasta --seqid foo --start 10 --end 20 > start.fa")

def test_17(capsys):
    run("bio convert ncov --fasta --type CDS -end 10 > cds.fa")

def test_18(capsys):
    run("bio convert ncov --fasta --type gene --end 10 > start-gene.fa")

def test_19(capsys):
    run("bio convert ncov --fasta --translate --type CDS --end 10 > translates.fa")

def test_20(capsys):
    run("bio convert ncov --fasta --protein --start -10 > protein-end.fa")

def test_21(capsys):
    run("bio convert ncov --fasta --type CDS --gene S --end 10 > cds-gene.fa")

def test_22(capsys):
    run("bio convert ncov:S --fasta --end 10 >  cds-gene.fa")

def test_23(capsys):
    run("bio convert ncov --id YP_009724390.1 --fasta --end 10 >  cds-gene.fa")

def test_24(capsys):
    run("bio convert ncov:S --fasta --protein --seqid foo > cds-prot.fa")

def test_25(capsys):
    run("bio convert ATGGGC -i --fasta > inter.fa")

def test_26(capsys):
    run("bio convert ATGGGC -i --translate --seqid foo >  inter-trans.fa")

def test_27(capsys):
    run("bio convert ATGGGC -i --revcomp --translate --seqid foo > inter-revcomp1.fa")

def test_28(capsys):
    run("bio convert ATGGGC -i --reverse --complement --translate --seqid foo >  inter-revcomp2.fa")

def test_29(capsys):
    run("bio align ncov ratg13 --end 180 > align-dna.txt")

def test_30(capsys):
    run("bio align ncov:S ratg13:S --end 180 > align-gene.txt")

def test_31(capsys):
    run("bio align ncov:S ratg13:S --end 180 --translate > align-translation.txt")

def test_32(capsys):
    run("bio align ratg13:S ncov:S  --end 180 -1 > align-gene-pept1.txt")

def test_33(capsys):
    run("bio align ratg13:S ncov:S  --end 180 -3 > align-gene-pept3.txt")

def test_34(capsys):
    run("bio align ncov:S ratg13:S --protein --end 60 > align-protein.txt")

def test_35(capsys):
    run("bio align THISLINE ISALIGNED  -i --local > align-local.txt")

def test_36(capsys):
    run("bio align THISLINE ISALIGNED -i --global > align-global.txt")

def test_37(capsys):
    run("bio align THISLINE ISALIGNED -i --semiglobal > align-semiglobal.txt")

def test_38(capsys):
    run("bio taxon 9606 > taxon_9606.txt")

def test_39(capsys):
    run("bio taxon 9606 --lineage > taxon_9606_lineage.txt")

def test_40(capsys):
    run("bio taxon 9606 --depth 1 > taxon_9606_depth_1.txt")

def test_41(capsys):
    run("bio taxon  ncov ratg13 > taxon_data.txt")

def test_42(capsys):
    run("echo 9606 | bio taxon > taxon_ids.txt")

def test_43(capsys):
    run("bio data --delete ebola")

def test_44(capsys):
    run("bio fetch KM233118 --rename ebola")

def test_45(capsys):
    run("bio runinfo ebola > sra-test.txt")

