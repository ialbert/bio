"""
This package attempts to simplify the BioPython SeqRecords into a simpler, flatter structure that
can be more readily worked with.
"""
import sys, os
from pprint import pprint
from intervaltree import Interval, IntervalTree
from biorun import utils

try:
    from Bio import SeqIO
except ImportError as exc:
    print(f"*** {exc}", file=sys.stderr)
    print(f"*** This software requires biopython.", file=sys.stderr)
    print(f"*** Try: conda install biopython", file=sys.stderr)
    sys.exit(1)


def get_attr(obj, name, default=''):
    """
    Returns an attribute stored in a BioPython record.
    """
    if isinstance(obj, dict):
        value = obj.get(name)
    elif hasattr(obj, name):
        value = getattr(obj, name)
    elif hasattr(obj, "qualifiers"):
        value = obj.qualifiers.get(name, default)
    elif hasattr(obj, "annotations"):
        value = obj.annotations.get(name, default)
    else:
        value = default

    return value


class Feature(object):
    """
    A thin, a more usable wrapper around a BioPython feature.

    Each feature has a type, has qualifiers and start, and stop coordinates.
    """

    def __init__(self, feat):
        """
        :param feat:  a Biopython feature object.
        """

        # The BioPython feature.
        self.feat = feat
        self.type = get_attr(self.feat, "type")

        self.qualifiers = dict(self.feat.qualifiers)

        # One based coordinate system.
        self.start = int(self.feat.location.start) + 1
        self.end = int(self.feat.location.end)

    def __repr__(self):
        """
        String representation of the feature.
        """
        return f"type:{self.type} loc:[{self.start}, {self.end}]"

    def __getattr__(self, item):
        """
        Attributes looked up in the qualifiers.
        """
        if item not in self.qualifiers:
            # print(self.qualifiers)
            raise AttributeError(item)

        return get_attr(self.feat, item)

    def as_dict(self):
        """
        Represents a Biopython Feature as a dictionary.
        """
        data = dict(
            start=self.start,
            end=self.end,
            type=self.type,
        )
        data.update(self.qualifiers)
        return data



class Sequence(object):
    """
    A simplified BioPython sequence record representation.

    Features from the original BioPython sequence record are stored as an interval tree of simplified features.
    """
    DEFAULT_ATTRS = {
        "id", "name", "dbxrefs", "description",
    }

    def __init__(self, rec):
        """
        A simpler, and more usable wrapper around a BioPython sequence record.
        """
        # Keep the Biopython record.
        self.rec = rec

        # Shortcut to annotations.
        self.ann = rec.annotations

        # Attributes that are always set on an instance.
        for attr in self.DEFAULT_ATTRS:
            value = get_attr(self.rec, attr)
            setattr(self, attr, value)

        # Store the features as intervals.
        self.features = IntervalTree()

        # Populate the interval tree.
        for feat in rec.features:
            obj = Feature(feat=feat)
            self.features[obj.start: obj.end] = obj

    def __iter__(self):
        """
        Iterates over the feature interval tree..
        """
        return iter(self.features)

    def __getattr__(self, item):
        """
        Attributes looks up in the annotations.
        """
        if item not in self.ann:
            raise AttributeError(item)

        return get_attr(self.ann, item)

    def as_dict(self):
        """
        Formats the sequence as dictionary data.
        """
        data = dict(
            annotations=dict(self.ann),
        )
        return data


def parse_genbank(stream, attrs=[]):
    stream = SeqIO.parse(stream, "genbank")
    recs = map(lambda rec: Sequence(rec), stream)
    return recs


def parse_fasta(stream, attrs=[]):
    recs = SeqIO.parse(stream, "fasta")
    return recs

def parse_file(fname, ftype=''):
    ftype = ftype or utils.guess_type(fname)

    if ftype == utils.GENBANK:
        recs = parse_genbank(fname)
    else:
        recs = parse_fasta(fname)

    return recs


if __name__ == "__main__":
    import doctest

    doctest.testmod()