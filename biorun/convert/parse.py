"""
Convert data from Entrez into different formats
"""
import os
import plac
from Bio import SeqIO
from biorun import utils
from biorun.data import fetch

EMBL, FASTA, GB = 'embl', 'fasta', 'genbank'


@utils.timer
def converter(in_file, in_format, out_file, out_format):
    SeqIO.convert(in_file=in_file,
                  in_format=in_format,
                  out_file=out_file,
                  out_format=out_format)
    return


@utils.timer
def fasta_converter(in_file, in_format, out_file, verb=0):
    """
    Convert a given file with a format to
    """
    # Load and parse the file.
    input_handle = open(in_file, 'r')
    stream = SeqIO.parse(input_handle, in_format)

    # Make sure the file is in fasta format.
    formatter = lambda seq_record: f">{seq_record.id} {seq_record.description}\n{seq_record.seq}\n"

    # Download fasta file into out file
    utils.download(stream, outname=out_file, verb=verb, formatter=formatter)

    return out_file


@plac.pos('acc', "accession number")
@plac.opt('db', "input target database ")
@plac.opt('mode', "input mode")
@plac.opt('output_dir', "output directory")
@plac.opt('format', "output data format")
@plac.flg('update', "overwrite existing data")
@plac.opt('verbosity', "verbosity level")
def run(acc, db='nuccore', in_format='gb', mode='text', output_dir=None, format=FASTA, update=False,
        verbosity=0):

    # Get a local copy of the genbank file or download it using Entrez.
    in_file = fetch.save_or_get(acc=acc, db=db, format=in_format,
                                mode=mode,
                                output_dir=output_dir,
                                update=update,
                                verbosity=verbosity)

    # Resolve the converted output name
    outname = utils.resolve_fname(acc=acc, directory=output_dir, ext=format)

    # Do not override existing file without update flag.
    if not update and os.path.isfile(outname):
        msg = f"File already exists at: {outname}"
        utils.print_message(msg=msg, styles=[utils.BOLD], verb=verbosity)
        return

    # Convert input file and write into outname.
    #converted = converter(in_file=in_file, in_format=in_format, out_file=outname, out_format=format)
    converted = fasta_converter(in_file=in_file, in_format=in_format, out_file=outname, verb=0)

    return
