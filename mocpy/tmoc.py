#!/usr/bin/env python
# -*- coding: utf-8 -*

"""
tmoc.py

functions to read/write and manipulate TimeMocs

"""

from __future__ import absolute_import, division, print_function, unicode_literals
from . import py23_compat

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from astropy.time import Time
from astropy.table import Table

from .abstract_moc import AbstractMoc

__author__ = "Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr"


class TimeMoc(AbstractMoc):
    DAY_MICRO_SEC = 86400000000.

    def __init__(self):
        AbstractMoc.__init__(self)

    @classmethod
    def from_table(cls, table, t_column, format='decimalyear'):
        """
        Create a tmoc instance from an astropy.table.Table object. The user has to specify the name of the time column
        to consider.

        :param table: the table containing all the sources we want to build the TimeMoc
        :param t_column: the name of the column to consider for building the TimeMoc
        :param format: the format of the times written in this column. See astropy.time.Time
        for the list of all available formats.
        :return: the tmoc instance resulting from the sources of the table
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
        Load a tmoc from a CSV/ascii file containing two columns (t_min, t_max)

        This method uses astropy.table.Table class and its read method to parse CSV/Ascii files.

        :param path: path leading to the CSV/ascii formatted file.
        This file must contains two columns whose names are t_min and t_max.
        :param kwargs: options that will be passed to the astropy.time.Time instances (e.g. format, scale options)
        :return: a newly created tmoc built from the given csv file
        """

        time_moc = TimeMoc()

        table = Table.read(path)
        try:
            t_min = Time(table['t_min'], **kwargs)
            t_max = Time(table['t_max'], **kwargs)
        except KeyError:
            raise KeyError("'t_min' or 't_max' are not part of the table columns : {0}".format(table.info()))

        add_all_time_intervals = np.vectorize(time_moc.add_time_interval)
        add_all_time_intervals(time_start=t_min,
                               time_end=t_max)
        return time_moc

    def add_time(self, time):
        """
        Add a single time observation to the current tmoc

        :param time: time observation (astropy.time.Time instance)
        :return:
        """
        if not isinstance(time, Time):
            raise TypeError("You must pass astropy.time.Time instances to the add_time_interval"
                            "method")
        time_us_start = int(time.jd * TimeMoc.DAY_MICRO_SEC)
        time_us_end = time_us_start + 1

        self._interval_set.add((time_us_start, time_us_end))

    def add_time_interval(self, time_start, time_end):
        """
        Add a time interval of observation to the current tmoc

        :param time_start: starting time of the observation (astropy.time.Time instance)
        :param time_end: ending time of the observation (astropy.time.Time instance)
        :return:
        """
        if not isinstance(time_start, Time) or not isinstance(time_end, Time):
            raise TypeError("You must pass astropy.time.Time instances to the add_time_interval"
                            "method")

        if time_start >= time_end:
            raise ValueError('time_start must be < compared to the time_end')

        time_us_start = int(time_start.jd * TimeMoc.DAY_MICRO_SEC)
        time_us_end = int(time_end.jd * TimeMoc.DAY_MICRO_SEC) + 1

        self._interval_set.add((time_us_start, time_us_end))

    @property
    def total_duration(self):
        """
        Return the total duration covered by the temporal moc in us
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
        Return a percentage of fill between the min and max time the moc is defined.
        A value near 0 shows a sparse temporal moc (i.e. the moc does not cover a lot
        of time and covers very distant times. A value near 1 means that the moc covers
        a lot of time without big pauses.
        """
        return self.total_duration / float(self.max_time - self.min_time)

    @property
    def min_time(self):
        """
        Return the min time of the tmoc in jd
        """
        return self._interval_set.min / TimeMoc.DAY_MICRO_SEC

    @property
    def max_time(self):
        """
        Return the max time of the tmoc in jd
        """
        return self._interval_set.max / TimeMoc.DAY_MICRO_SEC

    def filter_table(self, table, t_column, format='decimalyear', keep_inside=True):
        """
        Filter an astropy.table.Table to keep only rows inside (or outside) the tmoc instance

        :param table: the full astropy table to filter
        :param t_column: the name of the time column to consider for filtering the table
        :param format: the format of the time column
        :param keep_inside: boolean describing if we want to keep rows inside or outside the tmoc
        :return: the filtered Table (astropy.table)
        """

        time_arr = Time(table[t_column], format=format, scale='tai')
        pix_arr = (time_arr.jd * TimeMoc.DAY_MICRO_SEC)
        pix_arr = pix_arr.astype(int)

        intervals_arr = np.array(self._interval_set.intervals)
        inf_arr = np.vstack([pix_arr[i] >= intervals_arr[:, 0] for i in range(pix_arr.shape[0])])
        sup_arr = np.vstack([pix_arr[i] <= intervals_arr[:, 1] for i in range(pix_arr.shape[0])])

        if keep_inside:
            res = inf_arr & sup_arr
            filtered_rows = np.any(res, axis=1)
        else:
            res = ~inf_arr | ~sup_arr
            filtered_rows = np.all(res, axis=1)

        return table[filtered_rows]

    def add_fits_header(self, tbhdu):
        """
        Add TIMESYS to the header of the fits describing the tmoc

        :param tbhdu: fits python object
        """
        tbhdu.header['TIMESYS'] = ('JD', 'ref system JD BARYCENTRIC TT, 1 microsec level 29')

    def plot(self, title='TimeMoc', view=(None, None)):
        """
        Plot the tmoc with matplotlib

        :param title: the title of the plot
        :param view: a tuple (x_min, x_max) defining the time window of the plot
        """
        if self._interval_set.empty():
            print('Nothing to print. This TimeMoc object is empty.')
            return

        plot_order = 15
        if self.max_order > plot_order:
            plotted_moc = self.degrade_to_order(plot_order)
        else:
            plotted_moc = self

        plotted_moc._interval_set.intervals

        min_jd, max_jd = view

        if not min_jd:
            min_jd = plotted_moc.min_time
        if not max_jd:
            max_jd = plotted_moc.max_time

        if max_jd < min_jd:
            raise ValueError("Invalid selection: max_jd = {0} must be > to min_jd = {1}".format(max_jd, min_jd))

        fig1 = plt.figure(figsize=(9.5, 5))
        ax = fig1.add_subplot(111)

        ax.set_xlabel('iso')
        ax.get_yaxis().set_visible(False)

        size = 2000
        delta = (max_jd - min_jd) / size
        min_jd_time = min_jd

        ax.set_xticks([0, size])
        ax.set_xticklabels(Time([min_jd_time, max_jd], format='jd', scale='tdb').iso, rotation=70)

        y = np.zeros(size)
        for (s_time_us, e_time_us) in plotted_moc._interval_set.intervals:
            s_index = int((s_time_us / TimeMoc.DAY_MICRO_SEC - min_jd_time) / delta)
            e_index = int((e_time_us / TimeMoc.DAY_MICRO_SEC - min_jd_time) / delta)
            y[s_index:(e_index+1)] = 1.0

        # hack in case of full time mocs.
        if np.all(y):
            y[0] = 0

        z = np.tile(y, (int(size//10), 1))

        plt.title(title)

        color_map = LinearSegmentedColormap.from_list('w2r', ['#fffff0', '#aa0000'])
        color_map.set_under('w')
        color_map.set_bad('gray')

        plt.imshow(z, interpolation='bilinear', cmap=color_map)

        def on_mouse_motion(event):
            for txt in ax.texts:
                txt.set_visible(False)

            text = ax.text(0, 0, "", va="bottom", ha="left")

            time = Time(event.xdata * delta + min_jd_time, format='jd', scale='tdb')

            tx = '{0}'.format(time.iso)
            text.set_position((event.xdata - 50, 700))
            text.set_rotation(70)
            text.set_text(tx)

        cid = fig1.canvas.mpl_connect('motion_notify_event', on_mouse_motion)

        plt.show()
