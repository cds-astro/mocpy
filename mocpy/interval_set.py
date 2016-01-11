#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Manages a set of intervals
"""

from __future__ import print_function, division

# Python 3 support
try:
    xrange
except NameError:
    xrange = range

class IntervalSet:
    def __init__(self, interval_set=None):
        self.clear()
        if interval_set:
            self._intervals = list(interval_set._intervals)

    def clear(self):
        """
        Remove all entries
        """
        self._intervals = []
        self._must_check_consistency = False
        
    def empty(self):
        """
        Return True if the set is empty
        False otherwise
        """
        return len(self._intervals)==0
    
    def n_intervals(self):
        """
        Return the number of intervals in the set
        """
        return len(self._intervals)

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
            
            if itv[0] > stop: # gap between intervals
                ret.append((start, stop))
                start, stop = itv
            else: # merge intervals
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
        if self._must_check_consistency:
            self._intervals = IntervalSet.merge_intervals(self._intervals)
            self._must_check_consistency = False

        return self._intervals

    def add(self, item):
        """
        add a new item (integer or interval)
        """
        self._must_check_consistency = True
        if isinstance(item, tuple):
            if item[0]>item[1]:
                return
            self._intervals.append(item)

        else:
            self._intervals.append((item, item+1))
            
    def union2(self, another_is):
        """
        return the union between 2 IntervalSet
        """
        res = IntervalSet()
        if self.n_intervals()==0:
            res._intervals = another_is.intervals
        elif another_is.n_intervals()==0:
            res._intervals = self.intervals
        else:
            res._intervals =  IntervalSet.merge(self.intervals, another_is.intervals, lambda in_a, in_b: in_a or in_b)
        
        return res
    
    
    def union(self, another_is):
        """
        return the union between 2 IntervalSet
        """
        interval_set = IntervalSet()
        interval_set._intervals = IntervalSet.merge_intervals(self.intervals + another_is.intervals)
        
        return interval_set
    
    def difference(self, another_is):
        """
        return the difference between the current instance and another_is
        """
        res = IntervalSet()
        res._intervals =  IntervalSet.merge(self.intervals, another_is.intervals, lambda in_a, in_b: in_a and not in_b)
        
        return res
    
    def intersection(self,another_is):
        """
        return the intersection between the current instance and another_is
        """
        res = IntervalSet()
        res._intervals =  IntervalSet.merge(self.intervals, another_is.intervals, lambda in_a, in_b: in_a and in_b)
        
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
                 for i in xrange(0, len(list_of_endpoints) - 1, 2)]

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


        
        
        
            
if __name__=='__main__':
    iset = IntervalSet()
    import random
    for x in xrange(0, 10):
        start = random.randint(0, 100)
        iset.add((start, start+random.randint(0, 30)))

    another_iset = IntervalSet()
    for x in xrange(0, 10):
        start = random.randint(0, 100)
        another_iset.add((start, start+random.randint(0, 30)))

    print('A', iset.intervals)

    print('B', another_iset.intervals)
                     
    print('A union B', iset.union(another_iset).intervals)
    
    print('A - B', iset.difference(another_iset).intervals)
    
    print('B - A', another_iset.difference(iset).intervals)
    
    print('A intersection B', iset.intersection(another_iset).intervals)
    
    print(IntervalSet.flatten(iset.intervals))
    
    
