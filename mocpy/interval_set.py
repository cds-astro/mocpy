#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
interval_set.py

Manages a set of intervals

"""

from __future__ import absolute_import, division, print_function, unicode_literals
from .import py23_compat

__author__ = "Thomas Boch"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr"


class IntervalSet:
    def __init__(self, interval_set=None, intervals_l=None):
        self.clear()
        if interval_set:
            self._intervals = list(interval_set._intervals)
        if intervals_l:
            for item in intervals_l:
                self.add(item)

    def __repr__(self):
        return "{0}".format(self._intervals)

    def __eq__(self, another_is):
        if not isinstance(another_is, IntervalSet):
            raise TypeError
        return self._intervals == another_is._intervals

    @property
    def max(self):
        sorted_intervals = sorted(self._intervals, key=lambda interval: interval[1])
        return (sorted_intervals[-1])[1]

    @property
    def min(self):
        sorted_intervals = sorted(self._intervals, key=lambda interval: interval[0])
        return (sorted_intervals[0])[0]

    def clear(self):
        """
        Remove all entries
        """
        self._intervals = []
        self.__must_check_consistency = False

    def empty(self):
        """
        Return True if the set is empty
        False otherwise
        """
        return len(self._intervals) == 0

    @staticmethod
    def merge_intervals(intervals):
        """
        Merge adjacent intervals
        """
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

        return ret

    @property
    def intervals(self):
        """
        return the sorted list of intervals (merged if needed)
        """
        if self.__must_check_consistency:
            self._intervals = IntervalSet.merge_intervals(self._intervals)
            self.__must_check_consistency = False

        return self._intervals

    def add(self, item):
        """
        add a new item (integer or interval)
        """
        self.__must_check_consistency = True
        if isinstance(item, tuple):
            if item[0] > item[1]:
                return
            self._intervals.append(item)

        else:
            self._intervals.append((item, item+1))

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

    @staticmethod
    def flatten(interval_set):
        """Convert a list of intervals to a list of endpoints"""
        res = []
        for iv in interval_set:
            res.append(iv[0])
            res.append(iv[1])

        return res

    @staticmethod
    def unflatten(list_of_endpoints):
        """Convert a list of endpoints, with an optional terminating sentinel,
         into a list of intervals"""
        return [(list_of_endpoints[i], list_of_endpoints[i + 1])
                for i in range(0, len(list_of_endpoints) - 1, 2)]

    @staticmethod
    def merge(a_intervals, b_intervals, op):
        """Merge two lists of intervals according to the boolean function op"""
        a_endpoints = IntervalSet.flatten(a_intervals)
        b_endpoints = IntervalSet.flatten(b_intervals)

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

        return IntervalSet.unflatten(res)

