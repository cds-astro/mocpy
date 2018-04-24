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
    def from_table(cls, table, t_column, format='decimalyear'):
        """
        Create a TimeMoc from an astropy.table.Table object
        The user has to specify the column holding the time

        :param table: the table containing all the sources we want to build the TimeMoc
        :param t_column: the name of the column to consider for building the TimeMoc
        :param format: the format of the times written in this column
        :return: the TimeMoc object resulting from the sources of the table
        """
        if not isinstance(t_column, str):
            raise TypeError

        moc = TimeMoc()

        times = Time(table[t_column], format=format, scale='tai')

        for time_start in times:
            moc.add_time(time_start)

        return moc

    @classmethod
    def from_csv_file(cls, path, **kwargs):
        """
        Load a time moc from a CSV file containing two columns (t_min, t_max)

        :param path: path of the CSV file
        :param kwargs: args that will be passed to the astropy.time.Time objects (e.g. format, scale)
        :return: a time_moc
        """
        import csv

        time_moc = TimeMoc()

        with open(path, 'r') as csv_file:
            spam_reader = csv.reader(csv_file, delimiter=',', quotechar='|')
            row_num = 0
            for row in spam_reader:
                if row_num > 0:
                    t_min = Time(float(row[0]), **kwargs)
                    t_max = Time(float(row[1]), **kwargs)
                    time_moc.add_time_interval(t_min, t_max)

                row_num += 1

        return time_moc

    @classmethod
    def _from_specific_file(cls, hdulist, moc_order, path):
        # Entering this method triggers an error for TimeMoc because TimeMoc can only
        # be load from a set of nuniq intervals (as defined in the method from_file of AbstractMoc cls)
        raise FileNotFoundError('Error founding/opening file {0:s}'.format(path))

    def add_time(self, time):
        if not isinstance(time, Time):
            raise TypeError("You must pass astropy.time.Time instances to the add_time_interval"
                            "method")
        time_us_start = long(time.jd * TimeMoc.DAY_MICRO_SEC)
        time_us_end = time_us_start + 1

        self._interval_set.add((time_us_start, time_us_end))

    def add_time_interval(self, time_start, time_end):
        if not isinstance(time_start, Time) or not isinstance(time_end, Time):
            raise TypeError("You must pass astropy.time.Time instances to the add_time_interval"
                            "method")

        if time_start >= time_end:
            raise ValueError('time_start must be < compared to the time_end')

        time_us_start = long(time_start.jd * TimeMoc.DAY_MICRO_SEC)
        time_us_end = long(time_end.jd * TimeMoc.DAY_MICRO_SEC) + 1

        self._interval_set.add((time_us_start, time_us_end))

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

    def filter_table(self, table, t_column, format='decimalyear', keep_inside=True):
        """
        Filter an astropy.table.Table to keep only rows inside (or outside) the MOC instance
        Return the (newly created) filtered Table
        """
        time_arr = Time(table[t_column], format=format, scale='tai')
        pix_arr = (time_arr.jd * TimeMoc.DAY_MICRO_SEC)
        pix_arr = pix_arr.astype(int)

        itvs = np.array(self._interval_set.intervals)
        inf_arr = np.vstack([pix_arr[i] >= itvs[:, 0] for i in range(pix_arr.shape[0])])
        sup_arr = np.vstack([pix_arr[i] <= itvs[:, 1] for i in range(pix_arr.shape[0])])

        res = inf_arr & sup_arr
        if not keep_inside:
            res = np.logical_not(res)

        filtered_rows = np.any(res, axis=1)

        return table[filtered_rows]

    def add_fits_header(self, tbhdu):
        tbhdu.header['TIMESYS'] = ('JD', 'ref system JD BARYCENTRIC TT, 1 microsec level 29')

    def plot(self, title='TimeMoc', min_jd=None, max_jd=None):
        assert not self._interval_set.empty(), ValueError('Empty time moc instance')

        plot_order = 15
        if self.max_order > plot_order:
            plotted_moc = self.degrade_to_order(plot_order)
        else:
            plotted_moc = self

        plotted_moc._interval_set.intervals

        if not min_jd:
            min_jd = plotted_moc.min_time
        if not max_jd:
            max_jd = plotted_moc.max_time

        assert max_jd > min_jd, ValueError('max_jd must be > to min_jd')

        fig1 = plt.figure(figsize=(15, 20))
        ax = fig1.add_subplot(111)

        ax.set_xlabel('jd')
        ax.get_yaxis().set_visible(False)

        size = 2000
        delta = (max_jd - min_jd) / size
        min_jd_time = min_jd

        num_ticks = 5
        ax.set_xticks([x for x in range(0, size+1, int(size/num_ticks))])
        ax.set_xticklabels([str(int(x * delta + min_jd_time)) for x in range(0, size+1, int(size/num_ticks))])

        y = np.zeros(size)
        for (s_time_us, e_time_us) in plotted_moc._interval_set.intervals:
            s_index = int((s_time_us / TimeMoc.DAY_MICRO_SEC - min_jd_time) / delta)
            e_index = int((e_time_us / TimeMoc.DAY_MICRO_SEC - min_jd_time) / delta)
            y[s_index:(e_index+1)] = 1.0

        z = np.tile(y, (int(size//10), 1))

        plt.title(title)
        from matplotlib.colors import LinearSegmentedColormap
        color_map = LinearSegmentedColormap.from_list('w2r', ['#ffffff', '#aa0000'])
        color_map.set_under('w')
        color_map.set_bad('gray')

        plt.imshow(z, interpolation='bilinear', cmap=color_map)
        plt.show()

