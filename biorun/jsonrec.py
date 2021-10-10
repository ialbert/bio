import sys, json
from biorun.libs import placlib as plac

from biorun import parser

from biorun.parser import SOURCE, FEATURES, RECORD, ID, TYPE, ANNOTATIONS, LOCATIONS, SEQUENCE


def as_json(rec):
    data = {
        ID: rec.id, TYPE: rec.type, ANNOTATIONS: rec.ann, LOCATIONS: rec.locs
    }
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

            # Move the source to last
            if last:
                entry[FEATURES].append(last)

            # Create a new entry
            entry = {RECORD: rec.parent, FEATURES: [], SOURCE: str(rec.seq)}
            data.append(entry)
            last = as_json(rec)

        else:

            target = entry.get(FEATURES, [])
            row = as_json(rec)
            target.append(row)

    text = json.dumps(data, indent=4)

    print(text)
