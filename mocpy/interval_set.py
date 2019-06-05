# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
import copy
import numpy as np

from . import core

__author__ = "Thomas Boch"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr"


class IntervalSet:
    """Internal data structure for representing a MOC using the NESTED numbering scheme.

    A MOC object is a set of HEALPix cells at different orders (tuple (ipix, order)).
    MOC uses the NESTED numbering scheme and thus each HEALPix cell can be
    stored as one interval : [ipix*4^(29-order), (ipix+1)*4^(29-order)] - 29 being the maximum order of HEALPix cells
    one can encode in a 64 bit signed integer. See the `MOC IVOA standard paper <http://www.ivoa.net/documents/MOC/>`__
    for more explanations about how the NESTED numbering scheme works (i.e. how it is considered as a hierarchical
    numbering scheme).

    The key point is that a MOC can be represented by only one set of intervals (property of the NESTED scheme) that we
    are calling a consistent form i.e. such that it remains no overlapping intervals.
    An IntervalSet object must therefore keep its internal structure as a consistent set of intervals.
    For simplicity, the consistency step (i.e. the merge of the overlapping intervals) is done only once
    in the constructor. As there are no ways of modifying an IntervalSet object (e.g. add new HEALPix cells) then we are
    sure an IntervalSet is consistent when manipulating it for intersecting MOCs, doing their union etc...
    """
    HPY_MAX_ORDER = np.uint8(29)

    def __init__(self, intervals=None, make_consistent=True, min_depth=None):
        """
        IntervalSet constructor.

        The merging step of the overlapping intervals is done here.

        Parameters
        ----------
        intervals : `~numpy.ndarray`
            a N x 2 numpy array representing the set of intervals.
        make_consistent : bool, optional
            True by default. Remove the overlapping intervals that makes
            a valid MOC (i.e. can be plot, serialized, manipulated).
        """
        intervals = np.array([[]], dtype=np.uint64) if intervals is None else intervals
        # TODO: remove the cast to np.uint64
        # This code is executed as long as the Intervals objects
        # are not created from the rust code! (e.g. for the TimeMOCs)
        if intervals.dtype is not np.uint64:
            intervals = intervals.astype(np.uint64)
        self._intervals = intervals

        if make_consistent:
            if min_depth is None:
                min_depth = -1

            min_depth = np.int8(min_depth)
            self._merge_nested_intervals(min_depth)

    @classmethod
    def from_uniq(cls, pix):
        intervals = core.to_nested(pix)
        return cls(intervals=intervals, make_consistent=False)

    def _merge_nested_intervals(self, min_depth):
        if not self.empty():
            self._intervals = core.merge_nested_intervals(self._intervals, min_depth)

    def copy(self):
        """
        Deepcopy of self.

        Returns
        -------
        interval : `IntervalSet`
            a copy of self
        """
        return copy.deepcopy(self)

    def __repr__(self):
        return "{0}".format(self._intervals)

    def __eq__(self, another_is):
        """
        Equality operator override

        Parameters
        ----------
        another_is : `IntervalSet`
            IntervalSet object at the right of the equal operator

        Returns
        -------
        is_equal : bool
            boolean telling if self and ``another_is`` are equal or not.
        """
        return np.all(self._intervals == another_is._intervals)

    @property
    def min(self):
        return self._intervals.min()

    @property
    def max(self):
        return self._intervals.max()

    def empty(self):
        """
        Return True if the set is empty False otherwise
        """
        return self._intervals.size == 0

    def union(self, other):
        intervals = core.union(self._intervals, other._intervals)
        return IntervalSet(intervals, make_consistent=False)

    def intersection(self, other):
        intervals = core.intersection(self._intervals, other._intervals)
        return IntervalSet(intervals, make_consistent=False)

    def difference(self, other):
        intervals = core.difference(self._intervals, other._intervals)
        return IntervalSet(intervals, make_consistent=False)

    def complement(self):
        intervals = core.complement(self._intervals)
        return IntervalSet(intervals, make_consistent=False)

    @property
    def uniq(self):
        if self.empty():
            return self._intervals
        return core.to_uniq(self._intervals)

    @property
    def nested(self):
        return self._intervals