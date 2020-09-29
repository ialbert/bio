"""
This package attempts to simplify the BioPython SeqRecords into a simpler, flatter structure that
can be more readily worked with.
"""
import sys, os
from pprint import pprint
from intervaltree import Interval, IntervalTree
from biorun import utils
from functools import lru_cache
from itertools import *

try:
    from Bio import SeqIO
except ImportError as exc:
    print(f"*** {exc}", file=sys.stderr)
    print(f"*** This software requires biopython.", file=sys.stderr)
    print(f"*** Try: conda install biopython", file=sys.stderr)
    sys.exit(1)

logger = utils.logger

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

    Each feature has a type, has attributes, start, and stop coordinates.
    """

    def __init__(self, feat):
        """
        :param feat:  a Biopython feature object.
        """

        # The BioPython feature.
        self.feat = feat
        self.type = get_attr(self.feat, "type")
        self.strand = feat.strand
        self.attr = dict(self.feat.qualifiers)

        # One based coordinate system.
        self.start = int(self.feat.location.start)
        self.end = int(self.feat.location.end)

    @property
    @lru_cache()
    def name(self):
        """
        Attempts to produce a meaningful name for the feature.
        """
        value = self.get("gene") or self.type
        return value

    def get(self, name):
        if hasattr(self, name):
            return getattr(self, name)
        value = self.attr.get(name, '')
        return utils.flatten(value)

    def as_gff(self, anchor):
        """
        Return a 11 field GFF ready list
        """
        strand = "+" if self.strand else "-"

        pairs = [ ("Name", "name"), ("Function", "function"), ("ID", "protein_id")]
        pairs = [ (k, self.get(v)) for k,v in pairs ]
        pairs = filter(lambda x: x[1], pairs)
        attrib = map(lambda x: f"{x[0]}={x[1]}", pairs)
        attrib = ";".join(attrib)

        data = [
            anchor, ".", self.type, self.start+1, self.end, ".", strand, ".", attrib
        ]
        return data

    def as_dict(self):
        """
        Returns the feature as a dictionary
        """
        data = dict(name=self.name, start=self.start, end=self.end, type=self.type)
        data.update(self.attr)
        return data

    def __str__(self):
        return str(self.as_dict())


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
        self.tree = IntervalTree()

        # Populate the interval tree.
        for feat in rec.features:
            obj = Feature(feat=feat)
            self.tree[obj.start: obj.end] = obj

    def features(self, start=1, end=None, name=None, typ=None):
        """
        Returns the list of features.
        """

        # Slice the interval tree by range.
        logger.info("slicing the interval tree")
        feats = self.tree[start-1:end]

        # Sort by coordinates.
        logger.info("sorting features")
        feats = sorted(feats)

        # Return the features stored under the interval tree.
        feats = map(lambda f: f.data, feats)

        # Filter by name if set.
        feats = filter(lambda r: r.name == name, feats) if name else feats

        # Filter by type if requested.
        feats = filter(lambda r: r.type == typ, feats) if typ else feats

        return feats

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

    def __str__(self):
        return str(self.as_dict())


def parse_genbank(stream):
    recs = SeqIO.parse(stream, utils.GENBANK)
    recs = map(lambda rec: Sequence(rec), recs)

    return recs


if __name__ == "__main__":
    import doctest

    doctest.testmod()
