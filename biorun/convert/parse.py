"""
Convert data from Entrez into different formats
"""
import plac


def parse_fasta():
    """
    Given a file with the GenBank format, parse and return the fasta
    """
    return


def parse_gff():
    """
    Given a file with the GenBank format, parse and return the fasta
    """
    return


@plac.pos('acc', "accession number")
@plac.opt('db', "target database ")
@plac.opt('format', "output format")
@plac.opt('mode', "output mode")
@plac.opt('directory', "output directory to store data in")
@plac.opt('overwrite', "overwrite existing data when downloading.")
def run(acc, db='nuccore', format='gb', mode='text', directory=None, overwrite=False):

    # Check if the file has been downloaded locally first.
    return
    #stream = efetch(acc=acc, db=db, format=format, mode=mode)

    # Resolve file name from accession number and download.
    #outname = utils.resolve_fname(acc=acc, directory=directory)

    #utils.download(stream=stream, outname=outname, overwrite=overwrite)
