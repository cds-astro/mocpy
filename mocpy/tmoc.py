#!/usr/bin/env python
# -*- coding: utf-8 -*

"""

tmoc.py:
  functions to read/write and manipulate TimeMocs

"""

__author__ = "Baumann Matthieu"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

try:
    set
except NameError:
    from sets import Set as set
import numpy as np
from .interval_set import IntervalSet
from .moc import number_trailing_zeros

from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from astropy.time import Time, TimeDelta

import sys
if sys.version > '3':
    long = int

# Python 3 support
try:
    xrange
except NameError:
    xrange = range


class TimeMoc:
    VIZ_TABLE_MOC_ROOT_URL = ''
    VIZ_CAT_MOC_ROOT_URL = ''

    DAY_MICRO_SEC = 86400000000.
    HPY_MAX_NORDER = 29

    def __init__(self):
        self.__interval_set = IntervalSet()
        self._counter = 0

    @property
    def max_order(self):
        """
        This returns the deepest order needed to describe the current _interval_set
        """

        # TODO: cache value
        combo = self.total_duration

        ret = TimeMoc.HPY_MAX_NORDER - int(number_trailing_zeros(combo)//2)
        if ret < 0:
            ret = 0

        return ret

    def intersection(self, another_moc):
        """
        intersection with another MOC
        """
        iv_set_intersection = self.__interval_set.intersection(another_moc._interval_set)

        return TimeMoc.from_interval_set(iv_set_intersection)

    def union(self, another_moc):
        """
        union with another MOC
        """
        iv_set_union = self.__interval_set.union(another_moc._interval_set)

        return TimeMoc.from_interval_set(iv_set_union)

    @classmethod
    def from_interval_set(cls, interval_set):
        """
        Create a MOC from an IntervalSet (all intervals at deepest norder)
        """
        moc = TimeMoc()
        moc._interval_set = interval_set

        return moc

    def query_simbad(self, max_rows=10000):
        """
        query a view of SIMBAD data
        for SIMBAD objects in the coverage of the MOC instance
        """
        return self._query('SIMBAD', max_rows)

    def query_vizier_table(self, table_id, max_rows=10000):
        """
        query a VizieR table
        for sources in the coverage of the MOC instance
        """
        return self._query(table_id, max_rows)

    def _query(self, resource_id, max_rows):
        """
        internal method to query Simbad or a VizieR table
        for sources in the coverage of the MOC instance
        """
        pass

    def __ensure_consistency(self):
        """
        Force IntervalSet internal consistency
        """
        self.__interval_set.intervals

    def add_time_interval(self, time_start, time_end):
        if not isinstance(time_start, Time) or not isinstance(time_end, Time):
            raise TypeError("You must pass astropy.time.Time instances to the add_time_interval"
                            "method")

        if time_start >= time_end:
            raise ValueError('Starting time must be < to ending time')

        self._counter += 1
        # force consistency to prevent too large interval array
        if self._counter == 1000:
            self.__ensure_consistency()
            self._counter = 0

        self.__interval_set.add((time_start, time_end))

    @property
    def total_duration(self):
        """
        :return: the total duration covered by the temporal moc in microsec
        """
        if self.__interval_set.empty():
            return 0

        total_time = TimeDelta(0, format='jd')
        # The interval set is checked for consistency before looping over all the intervals
        for (start_time, stop_time) in self.__interval_set.intervals:
            total_time = total_time + (stop_time - start_time)

        total_time_us = total_time.jd * TimeMoc.DAY_MICRO_SEC

        if not total_time_us.is_integer():
            total_time_us += 1

        return int(total_time_us)

    @property
    def min_time(self):
        return self.__interval_set.min.jd

    @property
    def max_time(self):
        return self.__interval_set.max.jd

    def write(self, path, format='fits', optional_kw_dict=None):
        """
        Serialize a temporal moc in FITS in a given path
        Format can be 'fits' or 'json', though only the fits format is
        officially supported by the IVOA
        """
        formats = ('fits', 'json')
        if format not in formats:
            raise ValueError('format should be one of %s' % (str(formats)))

        # We ensure the temporal is consistent before dumping it in a file for save
        self.__ensure_consistency()

        if format == 'fits':
            moc_order = self.max_order
            if moc_order <= 13:
                format = '1J'
            else:
                format = '1K'

            """
            tbhdu = fits.BinTableHDU.from_columns(
                fits.ColDefs([fits.Column(name='UNIQ', format=format, array=np.array(uniq_array))]))
            tbhdu.header['PIXTYPE'] = 'HEALPIX'
            tbhdu.header['ORDERING'] = 'NUNIQ'
            tbhdu.header['COORDSYS'] = 'C'
            tbhdu.header['MOCORDER'] = moc_order
            tbhdu.header['MOCTOOL'] = 'MOCPy'
            if optional_kw_dict:
                for key, val in optional_kw_dict.items():
                    tbhdu.header[key] = val

            thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
            thdulist.writeto(path, clobber=True)
            """
        elif format == 'json':
            import json

            json_moc = {}
            total_time_us = self.total_duration

            min_order = TimeMoc.HPY_MAX_NORDER - TimeMoc.last_power_of_four(total_time_us)
            interval_set_tmp = IntervalSet(self.__interval_set)
            for order in range(min_order, self.max_order + 1):
                print(order)
                i_pix_l, res_interval = self.__ipix_in_interval_set(interval_set_tmp, order)
                if len(i_pix_l) > 0:
                    interval_set_tmp = interval_set_tmp.difference(res_interval)
                    json_moc[str(order)] = i_pix_l

            assert interval_set_tmp.empty(), ValueError('Error bad json construction : {0}'.format(json_moc))

            with open(path, 'w') as h:
                h.write(json.dumps(json_moc, sort_keys=True, indent=2))

    def __ipix_in_interval_set(self, interval_set, order):
        time_length_order_us = 1 << (2 * (TimeMoc.HPY_MAX_NORDER - order))
        max_time_us = int(self.max_time * TimeMoc.DAY_MICRO_SEC)
        min_time_us = int(self.min_time * TimeMoc.DAY_MICRO_SEC)

        num_time_len = min_time_us / time_length_order_us
        cur_time_us = int(num_time_len) * time_length_order_us
        if not num_time_len.is_integer():
            cur_time_us += time_length_order_us

        next_time_us = cur_time_us + time_length_order_us

        i_pix_l = []
        cur_i_pix = cur_time_us // time_length_order_us
        res_interval = IntervalSet()
        while cur_time_us < max_time_us:
            step = IntervalSet()
            step.add((Time(cur_time_us / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai'),
                      Time(next_time_us / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')))
            if interval_set.intersection(step) == step:
                i_pix_l.append(cur_i_pix)
                res_interval = res_interval.union(step)

            cur_i_pix += 1
            cur_time_us = next_time_us
            next_time_us += time_length_order_us

        return i_pix_l, res_interval

    @staticmethod
    def last_power_of_four(v):
        import math

        if v <= 0:
            raise ValueError('v must be > 0')

        if math.log(v, 2).is_integer():
            return int(math.log(v, 2) // 2)

        v -= 1
        v |= v >> 1
        v |= v >> 2
        v |= v >> 4
        v |= v >> 8
        v |= v >> 16
        v += 1

        power_of_two = math.log(v, 2) - 1
        return int(power_of_two // 2)

