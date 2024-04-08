import sys, json
from biorun.libs import placlib as plac
from pprint import pprint
from biorun import parser

from biorun.parser import SOURCE, FEATURES, RECORD, ID, TYPE, ANNOTATIONS
from biorun.parser import LOCATIONS, SEQUENCE, TITLE


def as_json(rec):
    data = {
        ID: rec.id, TYPE: rec.type,
    }

    data.update(rec.desc)

    data.update(
        {
            ANNOTATIONS: rec.ann, LOCATIONS: rec.locs
        }
    )

    return data


@plac.pos("fnames", "input files")
def run(*fnames):
    """
    Generates a json output.
    """
    # Get a stream of
    recs = parser.get_records(fnames)

    data = []

    entry = last = {}

    for rec in recs:

        if rec.type == SOURCE:
            # Create a new entry
            entry = {RECORD: rec.parent, FEATURES: [], SOURCE: str(rec.seq)}
            data.append(entry)

        target = entry.get(FEATURES, [])
        row = as_json(rec)
        target.append(row)

    text = json.dumps(data, indent=4)

    print(text)
