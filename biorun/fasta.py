"""
Handles FASTA related outputs
"""

import json
import sys

import plac
from biorun import jsonx
from biorun import utils

logger = utils.logger


@plac.flg("features", "convert the features")
@plac.opt("start", "start coordinate")
@plac.opt("end", "end coordinate")
@plac.opt("type_", "filter for a feature type")
@plac.opt("id_", "filter for a sequence id")
@plac.opt("name", "filter for a sequence name")
@plac.opt("gene", "filter for a gene name")
@plac.flg("proteins", "extract embedded protein sequences")
@plac.flg("translate", "translate DNA sequences", abbrev='R')
def run(features=False, proteins=False, translate=False,
        start='0', end='0', type_='', id_='', name='', gene=''):
    """
    Converts data to fastyaq
    """

    # Load the data
    data = json.loads(sys.stdin.read())

    # Parse start and end into user friendly numbers.
    start = utils.parse_number(start)
    end = utils.parse_number(end)
    ftype = type_
    seqid = id_

    recs = jsonx.select_records(data, features=features,
                                proteins=proteins, translate=translate,
                                start=start, end=end, ftype=ftype, seqid=seqid, name=name, gene=gene)

    for rec in recs:
        print(rec.format("fasta"), end='')
