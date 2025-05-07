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


@plac.flg("lines", "line-delimited JSON (default is a JSON array)", abbrev="l")
@plac.pos("fnames", "input files")
def run(lines=False, *fnames):
    """
    Generates a json output.
    """
    # Get a stream of
    recs = parser.get_records(fnames)

    data = []

    entry = last = {}

    for rec in recs:

        if rec.type == SOURCE:
            if lines and entry:
                print(json.dumps(entry, indent=None, separators=(',', ':')))

            # Create a new entry
            entry = {RECORD: rec.parent, FEATURES: [], SOURCE: str(rec.seq)}

            if not lines:
                data.append(entry)

        target = entry.get(FEATURES, [])
        row = as_json(rec)
        target.append(row)

    if lines:
        if entry:
            print(json.dumps(entry, indent=None, separators=(',', ':')))
    else:
        print(json.dumps(data, indent=4))
