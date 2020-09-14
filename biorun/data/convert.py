"""
Convert data from Entrez into different formats
"""
import os
import logging
import plac
from Bio import SeqIO
from biorun import utils
from biorun.data import fetch

from BCBio import GFF

logger = logging.getLogger('bio')

GFF3, EMBL, FASTA, GB = 'gff', 'embl', 'fasta', 'genbank'


def fasta_formatter(seq_record):
    return f">{seq_record.id} {seq_record.description}\n{seq_record.seq}\n"


def gff_formatter(seq_record):
    return seq_record


CONVERTERS = {FASTA: fasta_formatter,
              GFF3: gff_formatter}


def save_gff(sequence_stream, out_file):

    out_handle = open(out_file, "w")

    GFF.write(sequence_stream, out_handle)

    out_handle.close()

    return out_file


def converter(in_file, in_format, out_file, out_format=FASTA):
    """
    Convert a given file with a format to
    """
    # Load and parse the file.
    input_handle = open(in_file, 'r')

    sequence_stream = SeqIO.parse(input_handle, in_format)

    # Download stream into out file
    if out_format == GFF3:
        save_gff(in_file, in_format)
    else:
        # Fetch the formatter.
        formatter = CONVERTERS[out_format]
        utils.save_file(sequence_stream, outname=out_file, formatter=formatter)

    return out_file


@plac.pos('acc', "accession number")
@plac.opt('db', "input target database ")
@plac.opt('mode', "input mode")
@plac.opt('dir', "output directory", abbrev='output_dir')
@plac.opt('name', "output file name ( extension included )")
@plac.opt('format', "output data format")
@plac.flg('update', "overwrite existing data")
def run(acc, db='nuccore', in_format='gb', mode='text', output_dir=None, name=None, format=FASTA,
        update=False):

    # Get a local copy of the genbank file or save_file it using Entrez.
    in_file = fetch.save_or_get(acc=acc, db=db, format=in_format, mode=mode,
                                output_dir=output_dir, output_name=name,
                                update=update)

    # Resolve the converted output name
    outname = utils.resolve_fname(acc=acc, directory=output_dir, ext=format, output_name=name)

    # Do not override existing file without update flag.
    if not update and os.path.isfile(outname):
        logger.info(f"File already exists at: {outname}")
        return

    # Convert input file and write into outname.
    converter(in_file=in_file, in_format=in_format, out_file=outname, out_format=format)

    return
