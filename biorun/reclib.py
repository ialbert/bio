"""
Utility function for transforming BioPython SeqRecords
"""
import sys

from biorun import utils

logger = utils.logger

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError as exc:
    logger.error(f"{__name__}", stop=False)
    logger.error(f"{exc}", stop=False)
    logger.error(f"This software requires biopython.")
    logger.error(f"Try: conda install biopython>=1.78")
    sys.exit()


def filter_gene(rec, gene=None):
    """
    Filters BioPython SeqRecords by name
    """
    return True if not gene else rec.__gene__ == gene



def filter_name(rec, name=None):
    """
    Filters BioPython SeqRecords by name
    """
    #print (rec.name)
    return True if not name else (rec.name == name or rec.id == name)


def filter_seqid(rec, seqid=None):
    """
    Filters BioPython SeqRecords by id
    """
    return True if not seqid else rec.id == seqid


def filter_type(rec, ftype=None):
    """
    Filters BioPython SeqRecords by type.
    """
    return True if not ftype else rec.__biotype__ == ftype


def filter_id(rec, value=None):
    """
    Filters BioPython SeqRecords by id.
    """
    return True if not value else rec.id == id


def translate(rec):
    """
    Translates a BioPython Sequence record
    """
    # Truncate to nearest codon
    end = len(rec.seq) // 3 * 3
    rec.seq = rec.seq[:end].translate()
    rec.description = f"translated {rec.description}"
    return rec


def slice(rec, start=0, end=0):
    """
    Slices a BioPython Sequence record
    """

    # Nothing needs to be done
    if not (start or end):
        return rec

    # There is end but no start
    if end and not start:
        start = 1

    # There start but no end
    if start and not end:
        end = len(rec.seq)

    rec.seq = rec.seq[start - 1:end]
    rec.description = f"{start}:{end} {rec.description}"
    return rec

