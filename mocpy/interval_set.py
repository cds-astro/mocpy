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
        self.__must_check_consistency = True

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

    """
    @property
    def max(self):
        sorted_intervals = sorted(self._intervals, key=lambda interval: interval[1])
        return (sorted_intervals[-1])[1]

    @property
    def min(self):
        sorted_intervals = sorted(self._intervals, key=lambda interval: interval[0])
        return (sorted_intervals[0])[0]
    """
    def clear(self):
        """
        Remove all entries
        """
        self._intervals = np.array([])
        self.__must_check_consistency = False

    def empty(self):
        """
        Return True if the set is empty
        False otherwise
        """
        return self._intervals.size == 0

    @staticmethod
    def merge_intervals(intervals_arr):
        """
        Merge adjacent intervals
        """
        intervals = intervals_arr.tolist()

        ret = []
        start = stop = None
        for itv in sorted(intervals):
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

        return np.asarray(ret)

    @property
    def intervals(self):
        """
        return the sorted list of intervals (merged if needed)
        """
        if self.__must_check_consistency:
            self._intervals = IntervalSet.merge_intervals(self._intervals)
            self.__must_check_consistency = False

        return self._intervals

    def add_numpy_arr(self, np_arr):
        self.__must_check_consistency = True
        if self.empty():
            self._intervals = np_arr
            return
        self._intervals = np.vstack((self._intervals, np_arr))

    def add(self, item):
        """
        add a new item (i.e. interval)
        """
        self.__must_check_consistency = True
        if item[0] > item[1]:
            return
        if self.empty():
            self._intervals = np.array(item)
            return
        self._intervals = np.vstack((self._intervals, item))

    def union(self, another_is):
        """
        return the union between 2 IntervalSet
        """
        res = IntervalSet()
        if self.empty():
            res._intervals = another_is.intervals
        elif another_is.empty():
            res._intervals = self.intervals
        else:
            res._intervals = IntervalSet.merge(self.intervals, another_is.intervals, lambda in_a, in_b: in_a or in_b)

        return res

    def difference(self, another_is):
        """
        return the difference between the current instance and another_is
        """
        res = IntervalSet()
        if another_is.empty() or self.empty():
            res._intervals = self.intervals
        else:
            res._intervals = IntervalSet.merge(self.intervals, another_is.intervals, lambda in_a, in_b: in_a and not in_b)

        return res

    def intersection(self, another_is):
        """
        return the intersection between the current instance and another_is
        """
        res = IntervalSet()
        if not another_is.empty() and not self.empty():
            res._intervals = IntervalSet.merge(self.intervals, another_is.intervals, lambda in_a, in_b: in_a and in_b)

        return res

    @classmethod
    def to_nuniq_interval_set(cls, orderipix_itv):
        HPY_MAX_NORDER = 29

        r2 = orderipix_itv.copy()
        r3 = IntervalSet()
        res = IntervalSet()

        if r2.empty():
            return res

        order = 0

        while not r2.empty():
            shift = int(2 * (HPY_MAX_NORDER - order))
            ofs = (int(1) << shift) - 1
            ofs2 = int(1) << (2 * order + 2)

            r3.clear()
            for iv in r2.intervals:
                a = (int(iv[0]) + ofs) >> shift
                b = int(iv[1]) >> shift
                r3.add((a << shift, b << shift))
                res.add((a + ofs2, b + ofs2))

            if not r3.empty():
                r2 = r2.difference(r3)

            order += 1

        return res

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
        rtmp = IntervalSet()
        last_order = 0
        HPY_MAX_NORDER = 29
        intervals = nuniq_itv.intervals
        diff_order = HPY_MAX_NORDER
        shift_order = 2 * diff_order
        for interval in intervals:
            for j in range(interval[0], interval[1]):
                order, i_pix = uniq2orderipix(j)

                if order != last_order:
                    orderipix_itv = orderipix_itv.union(rtmp)
                    rtmp.clear()
                    last_order = order
                    diff_order = HPY_MAX_NORDER - order
                    shift_order = 2 * diff_order

                rtmp.add((i_pix << shift_order, (i_pix + 1) << shift_order))

        orderipix_itv = orderipix_itv.union(rtmp)
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

