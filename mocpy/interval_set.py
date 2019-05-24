# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
import copy
import numpy as np

from .py23_compat import range
from .utils import uniq2orderipix

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
        self._intervals = np.array([], dtype=np.uint64) if intervals is None else intervals

        if make_consistent:
            if min_depth is not None:
                min_depth = np.uint8(min_depth)
            self._merge_intervals(min_depth=min_depth)

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

    def _merge_intervals(self, min_depth):
        """
        Merge overlapping intervals.

        This method is called only once in the constructor.
        """
        def add_interval(ret, start, stop):
            if min_depth is not None:
                shift = np.uint8(2) * (np.uint8(29) - min_depth)
                mask = (np.uint64(1) << shift) - np.uint64(1)

                if stop - start < mask:
                    ret.append((start, stop))
                else:
                    ofs = start & mask
                    st = start
                    if ofs > 0:
                        st = (start - ofs) + (mask + np.uint64(1))
                        ret.append((start, st))

                    while st + mask + 1 < stop:
                        ret.append((st, st + mask + np.uint64(1)))
                        st = st + mask + np.uint64(1)

                    ret.append((st, stop))
            else:
                ret.append((start, stop))

        ret = []
        start = stop = None
        # Use numpy sort method
        self._intervals.sort(axis = 0)
        for itv in self._intervals:
            if start is None:
                start, stop = itv
                continue

            # gap between intervals
            if itv[0] > stop:
                add_interval(ret, start, stop)
                start, stop = itv
            else:
                # merge intervals
                if itv[1] > stop:
                    stop = itv[1]

        if start is not None and stop is not None:
            add_interval(ret, start, stop)

        self._intervals = np.asarray(ret)
        #self._intervals = self._intervals.astype(np.uint64)

    def union(self, another_is):
        """
        Return the union between self and ``another_is``.

        Parameters
        ----------
        another_is : `IntervalSet`
            an IntervalSet object.
        Returns
        -------
        interval : `IntervalSet`
            the union of self with ``another_is``.
        """
        result = IntervalSet()
        if another_is.empty():
            result._intervals = self._intervals
        elif self.empty():
            result._intervals = another_is._intervals
        else:
            # res has no overlapping intervals
            result._intervals = IntervalSet.merge(self._intervals,
                                                  another_is._intervals,
                                                  lambda in_a, in_b: in_a or in_b)
        return result

    def difference(self, another_is):
        """
        Return the difference between self and ``another_is``.

        Parameters
        ----------
        another_is : `IntervalSet`
            an IntervalSet object.
        Returns
        -------
        interval : `IntervalSet`
            the difference of self with ``another_is``.
        """
        result = IntervalSet()
        if another_is.empty() or self.empty():
            result._intervals = self._intervals
        else:
            result._intervals = IntervalSet.merge(self._intervals,
                                                  another_is._intervals,
                                                  lambda in_a, in_b: in_a and not in_b)
        return result

    def intersection(self, another_is):
        """
        Return the intersection between self and ``another_is``.

        Parameters
        ----------
        another_is : `IntervalSet`
            an IntervalSet object.
        Returns
        -------
        interval : `IntervalSet`
            the intersection of self with ``another_is``.
        """
        result = IntervalSet()
        if not another_is.empty() and not self.empty():
            result._intervals = IntervalSet.merge(self._intervals,
                                                  another_is._intervals,
                                                  lambda in_a, in_b: in_a and in_b)
        return result

    @classmethod
    def to_nuniq_interval_set(cls, nested_is):
        """
        Convert an IntervalSet using the NESTED numbering scheme to an IntervalSet containing UNIQ numbers for HEALPix
        cells.

        Parameters
        ----------
        nested_is : `IntervalSet`
            IntervalSet object storing HEALPix cells as [ipix*4^(29-order), (ipix+1)*4^(29-order)[ intervals.

        Returns
        -------
        interval : `IntervalSet`
            IntervalSet object storing HEALPix cells as [ipix + 4*4^(order), ipix+1 + 4*4^(order)[ intervals.
        """
        r2 = nested_is.copy()
        res = []

        if r2.empty():
            return IntervalSet()

        order = np.uint8(0)
        while not r2.empty():
            shift = np.uint8(2) * (np.uint8(IntervalSet.HPY_MAX_ORDER) - order)
            ofs = (np.uint64(1) << shift) - np.uint64(1)
            ofs2 = np.uint64(1) << (np.uint8(2) * order + np.uint8(2))

            r4 = []
            for iv in r2._intervals:
                a = (iv[0] + ofs) >> shift
                b = iv[1] >> shift

                c = a << shift
                d = b << shift
                if d > c:
                    r4.append((c, d))
                    res.append((a + ofs2, b + ofs2))

            if len(r4) > 0:
                r4_is = IntervalSet(np.asarray(r4, dtype=np.uint64))
                r2 = r2.difference(r4_is)

            order += np.uint8(1)

        return IntervalSet(np.asarray(res))

    @classmethod
    def from_nuniq_interval_set(cls, nuniq_is):
        """
        Convert an IntervalSet containing NUNIQ intervals to an IntervalSet representing HEALPix
        cells following the NESTED numbering scheme.

        Parameters
        ----------
        nuniq_is : `IntervalSet`
            IntervalSet object storing HEALPix cells as [ipix + 4*4^(order), ipix+1 + 4*4^(order)[ intervals.

        Returns
        -------
        interval : `IntervalSet`
            IntervalSet object storing HEALPix cells as [ipix*4^(29-order), (ipix+1)*4^(29-order)[ intervals.
        """
        nested_is = IntervalSet()
        # Appending a list is faster than appending a numpy array
        # For these algorithms we append a list and create the interval set from the finished list
        rtmp = []
        last_order = np.uint8(0)
        intervals = nuniq_is._intervals
        diff_order = np.uint8(IntervalSet.HPY_MAX_ORDER)
        shift_order = np.uint8(2) * diff_order
        for interval in intervals:
            for uniq in np.arange(np.int(interval[0]), np.int(interval[1]), dtype=np.uint64):
                order, i_pix = uniq2orderipix(uniq)

                if order != last_order:
                    nested_is = nested_is.union(IntervalSet(np.asarray(rtmp, dtype=np.uint64)))
                    rtmp = []
                    last_order = order
                    diff_order = np.uint8(IntervalSet.HPY_MAX_ORDER) - order
                    shift_order = np.uint8(2) * diff_order

                rtmp.append((i_pix << shift_order, (i_pix + np.uint64(1)) << shift_order))

        nested_is = nested_is.union(IntervalSet(np.asarray(rtmp, dtype=np.uint64)))
        return nested_is

    @staticmethod
    def merge(a_intervals, b_intervals, op):
        """
        Merge two lists of intervals according to the boolean function op

        ``a_intervals`` and ``b_intervals`` need to be sorted and consistent (no overlapping intervals).
        This operation keeps the resulting interval set consistent.

        Parameters
        ----------
        a_intervals : `~numpy.ndarray`
            A sorted merged list of intervals represented as a N x 2 numpy array
        b_intervals : `~numpy.ndarray`
            A sorted merged list of intervals represented as a N x 2 numpy array
        op : `function`
            Lambda function taking two params and returning the result of the operation between
            these two params.
            Exemple : lambda in_a, in_b: in_a and in_b describes the intersection of ``a_intervals`` and
            ``b_intervals`` whereas lambda in_a, in_b: in_a or in_b describes the union of ``a_intervals`` and
            ``b_intervals``.

        Returns
        -------
        array : `numpy.ndarray`
            a N x 2 numpy containing intervals resulting from the op between ``a_intervals`` and ``b_intervals``.
        """
        a_endpoints = a_intervals.flatten().tolist()
        b_endpoints = b_intervals.flatten().tolist()

        sentinel = max(a_endpoints[-1], b_endpoints[-1]) + 1

        a_endpoints += [sentinel]
        b_endpoints += [sentinel]

        a_index = 0
        b_index = 0

        res = []

        scan = min(a_endpoints[0], b_endpoints[0])
        while scan < sentinel:
            in_a = not ((scan < a_endpoints[a_index]) ^ (a_index % 2))
            in_b = not ((scan < b_endpoints[b_index]) ^ (b_index % 2))
            in_res = op(in_a, in_b)

            if in_res ^ (len(res) % 2):
                res += [scan]
            if scan == a_endpoints[a_index]:
                a_index += 1
            if scan == b_endpoints[b_index]:
                b_index += 1

            scan = min(a_endpoints[a_index], b_endpoints[b_index])

        return np.asarray(res, dtype=np.uint64).reshape((-1, 2))
