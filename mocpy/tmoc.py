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

import sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.time import Time
from astropy.io import fits

from .abstract_moc import AbstractMoc

if sys.version > '3':
    long = int

# Python 3 support
try:
    xrange
except NameError:
    xrange = range


class TimeMoc(AbstractMoc):
    DAY_MICRO_SEC = 86400000000.

    def __init__(self):
        AbstractMoc.__init__(self)

    @classmethod
    def _from_specific_file(cls, hdulist, moc_order, path):
        # Entering this method triggers an error for TimeMoc because TimeMoc can only
        # be load from a set of nuniq intervals (as defined in the method from_file of AbstractMoc cls)
        raise FileNotFoundError('Error founding/opening file {0:s}'.format(path))

    def add_time_interval(self, time_start, time_end):
        if not isinstance(time_start, Time) or not isinstance(time_end, Time):
            raise TypeError("You must pass astropy.time.Time instances to the add_time_interval"
                            "method")

        if time_start >= time_end:
            raise ValueError('time_start must be < compared to the time_end')

        """
        self.__counter += 1
        # force consistency to prevent too large interval array
        if self.__counter == 1000:
            self._ensure_consistency()
            self.__counter = 0
        """

        time_us_start = long(time_start.jd * TimeMoc.DAY_MICRO_SEC)
        time_us_end = time_end.jd * TimeMoc.DAY_MICRO_SEC
        if not time_us_end.is_integer():
            time_us_end += 1

        self._interval_set.add((time_us_start,
                                long(time_us_end)))

    @property
    def total_duration(self):
        """
        :return: the total duration covered by the temporal moc in us
        """
        if self._interval_set.empty():
            return 0

        total_time_us = 0
        # The interval set is checked for consistency before looping over all the intervals
        for (start_time, stop_time) in self._interval_set.intervals:
            total_time_us = total_time_us + (stop_time - start_time)

        return total_time_us

    @property
    def consistency(self):
        """
        :return: a percentage of fill between the min and max time the moc is defined.
        A value near 0 shows a sparse temporal moc (i.e. the moc does not cover a lot
        of time and covers very distant times. A value near 1 means that the moc covers
        a lot of time without big pauses.
        """
        return self.total_duration / float(self.max_time - self.min_time)

    @property
    def min_time(self):
        """Get the min time of the temporal moc in jd"""
        return self._interval_set.min / TimeMoc.DAY_MICRO_SEC

    @property
    def max_time(self):
        """Get the min time of the temporal moc in jd"""
        return self._interval_set.max / TimeMoc.DAY_MICRO_SEC

    def filter_table(self, table, keep_inside=True, format='decimalyear', *args):
        return self._filter(table,
                            keep_inside,
                            format,
                            *args)

    def _get_pix(self, row_values_l, n_side, format=None):
        assert len(row_values_l) == 1, ValueError('Filtering a table by a temporal moc is done on time'
                                                  ' columns such as u, g, r.. date')
        time = Time(row_values_l[0], format=format, scale='tai')

        i_pix = time.jd * TimeMoc.DAY_MICRO_SEC
        if not i_pix.is_integer():
            i_pix += 1

        return int(i_pix)

    def add_fits_header(self, tbhdu):
        tbhdu.header['TIMESYS'] = ('JD', 'ref system JD BARYCENTRIC TT, 1 microsec level 29')

    def plot(self):
        plot_order = 15
        if self.max_order > plot_order:
            plotted_moc = self.degrade_to_order(plot_order)
        else:
            plotted_moc = self

        plotted_moc._interval_set.intervals
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)

        ax.set_xlabel('jd')
        ax.get_yaxis().set_visible(False)

        size = 1000
        delta = (plotted_moc.max_time - plotted_moc.min_time) / size
        min_jd_time = plotted_moc.min_time

        num_ticks = 5
        ax.set_xticks([x for x in range(0, size+1, int(size/num_ticks))])
        ax.set_xticklabels([str(int(x * delta + min_jd_time)) for x in range(0, size+1, int(size/num_ticks))])

        y = np.zeros(size)
        for (s_time_us, e_time_us) in plotted_moc._interval_set.intervals:
            s_index = int((s_time_us / TimeMoc.DAY_MICRO_SEC - min_jd_time) / delta)
            e_index = int((e_time_us / TimeMoc.DAY_MICRO_SEC - min_jd_time) / delta)
            y[s_index:(e_index+1)] = 1.0

        z = np.tile(y, (int(size//10), 1))

        plt.imshow(z, interpolation='nearest')
        plt.show()

