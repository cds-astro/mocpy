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
    """Manages a set of intervals.

    TODO: describe this data structure,
    especially how things are stored.
    """
    def __init__(self, intervals=None):
        intervals = np.array([]) if intervals is None else intervals
        self._intervals = intervals
        self._merge_intervals()

    # TODO: from a Nx2 numpy array. To be removed lated
    @classmethod
    def from_numpy_array(cls, arr):
        return cls(arr)

    def copy(self):
        """Make a copy."""
        return copy.deepcopy(self)

    def __repr__(self):
        return "{0}".format(self._intervals)

    def __eq__(self, another_is):
        if not isinstance(another_is, IntervalSet):
            raise TypeError
        return np.all(self._intervals == another_is._intervals)

    def empty(self):
        """
        Return True if the set is empty
        False otherwise
        """
        return self._intervals.size == 0

    def _merge_intervals(self):
        """
        Merge adjacent intervals
        """

        ret = []
        start = stop = None
        for itv in sorted(self._intervals.tolist()):
            if start is None:
                start, stop = itv
                continue

            #  gap between intervals
            if itv[0] > stop:
                ret.append((start, stop))
                start, stop = itv
            else:
                #  merge intervals
                if itv[1] > stop:
                    stop = itv[1]

        if start is not None and stop is not None:
            ret.append((start, stop))

        self._intervals = np.asarray(ret)

    """@property
    def intervals(self):
        '''
        return the sorted list of intervals (merged if needed)
        '''
        if self.__must_check_consistency:
            self._intervals = IntervalSet.merge_intervals(self._intervals)
            self.__must_check_consistency = False

        return self._intervals
    """

    def union(self, another_is):
        """
        return the union between 2 IntervalSet
        """
        res = IntervalSet()
        if self.empty():
            another_is._merge_intervals()
            res._intervals = another_is._intervals
        elif another_is.empty():
            self._merge_intervals()
            res._intervals = self._intervals
        else:
            self._merge_intervals()
            another_is._merge_intervals()
            res._intervals = IntervalSet.merge(self._intervals, another_is._intervals, lambda in_a, in_b: in_a or in_b)

        # res has no overlapping intervals
        return res

    def difference(self, another_is):
        """
        return the difference between the current instance and another_is
        """
        res = IntervalSet()
        if another_is.empty() or self.empty():
            self._merge_intervals()
            res._intervals = self._intervals
        else:
            self._merge_intervals()
            another_is._merge_intervals()
            res._intervals = IntervalSet.merge(self._intervals, another_is._intervals, lambda in_a, in_b: in_a and not in_b)

        return res

    def intersection(self, another_is):
        """
        return the intersection between the current instance and another_is
        """
        res = IntervalSet()
        if not another_is.empty() and not self.empty():
            self._merge_intervals()
            another_is._merge_intervals()
            res._intervals = IntervalSet.merge(self._intervals, another_is._intervals, lambda in_a, in_b: in_a and in_b)

        return res

    @classmethod
    def to_nuniq_interval_set(cls, orderipix_itv):
        HPY_MAX_NORDER = 29

        r2 = orderipix_itv.copy()
        res = []

        if r2.empty():
            return res

        order = 0
        while not r2.empty():
            shift = int(2 * (HPY_MAX_NORDER - order))
            ofs = (int(1) << shift) - 1
            ofs2 = int(1) << (2 * order + 2)

            r4 = []
            for iv in r2._intervals:
                a = (int(iv[0]) + ofs) >> shift
                b = int(iv[1]) >> shift

                c = a << shift
                d = b << shift
                if d > c:
                    r4.append((c, d))

                res.append((a + ofs2, b + ofs2))

            if len(r4) > 0:
                r4_is = IntervalSet.from_numpy_array(np.asarray(r4))
                r2 = r2.difference(r4_is)

            order += 1

        return IntervalSet.from_numpy_array(np.asarray(res))

    @classmethod
    def from_nuniq_interval_set(cls, nuniq_itv):
        """
        Convert a nuniq interval set to a order/ipix interval set

        Parameters
        ----------
        nuniq_itv: `mocpy.interval_set.IntervalSet`
            a interval set object containing nuniq numbers

        Returns
        -------
        orderipix_itv: `mocpy.interval_set.IntervalSet`
            An order/ipix interval set
        """
        orderipix_itv = IntervalSet()
        # Appending a list is faster than appending a numpy array
        # For these algorithms we append a list and create the interval set from the finished list
        rtmp = []
        last_order = 0
        HPY_MAX_NORDER = 29
        intervals = nuniq_itv._intervals
        diff_order = HPY_MAX_NORDER
        shift_order = 2 * diff_order
        for interval in intervals:
            for j in range(interval[0], interval[1]):
                order, i_pix = uniq2orderipix(j)

                if order != last_order:
                    orderipix_itv = orderipix_itv.union(IntervalSet.from_numpy_array(np.asarray(rtmp)))
                    rtmp = []
                    last_order = order
                    diff_order = HPY_MAX_NORDER - order
                    shift_order = 2 * diff_order

                rtmp.append((i_pix << shift_order, (i_pix + 1) << shift_order))

        orderipix_itv = orderipix_itv.union(IntervalSet.from_numpy_array(np.asarray(rtmp)))
        return orderipix_itv

    @staticmethod
    def merge(a_intervals, b_intervals, op):
        """Merge two lists of intervals according to the boolean function op"""
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

        return np.asarray(res).reshape((-1, 2))
