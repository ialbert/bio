import plac

from Bio import SeqIO
from django.core.management import call_command
import logging
from biorun import utils
from biorun.data import fetch
from biorun.db.models import Locus,Feature

logger = logging.getLogger('bio')


def insert(input_file):

    input_handle = open(input_file, 'r')

    # Parse a stream of genbank file to be more of al
    stream = SeqIO.parse(input_handle, "gbk")

    for req in stream:
        print(dir(req))
        1/0

    # Iterate through the stream and insert into database.

    return


def migrate_db():
    """
    Make any necessary migrations and apply to database.
    """
    # Make new migrations first.
    call_command('makemigrations', 'db')

    # Apply new migrations.
    call_command('migrate', 'db')

    return


# Can also give a list of accession numbers or lists.
@plac.opt('acc', "comma separated list of accession numbers to add to database.")
@plac.flg('migrate', "migrate the database")
@plac.opt('db', "database to populate.")
def run(acc=None, migrate=False, db=None):
    """
    """

    if migrate:
        migrate_db()

    # Get the input files from using accession number
    # fname = fetch.get(acc=acc, db=db, format="gbk",
    #                 mode="text",
    #                 output_dir=input_dir)

    # Insert file into database.
    #insert(input_file=fname)

    return
