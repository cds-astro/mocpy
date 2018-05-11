#!/usr/bin/env python
# -*- coding: utf-8 -*

"""
tmoc.py

functions to read/write and manipulate TimeMocs

"""

from __future__ import absolute_import, division, print_function, unicode_literals
from . import py23_compat

import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from astropy.time import Time, TimeDelta
from astropy.table import Table

from .abstract_moc import AbstractMoc

__author__ = "Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr"


class TimeMoc(AbstractMoc):
    DAY_MICRO_SEC = 86400000000.
    # default observation time : 30 min
    DEFAULT_OBSERVATION_TIME = TimeDelta(30 * 60, format='sec', scale='tdb')

    def __init__(self):
        AbstractMoc.__init__(self)

    @classmethod
    def from_table(cls, table, t_column, delta_t=DEFAULT_OBSERVATION_TIME, **kwargs):
        """
        Create a TimeMoc instance from an `~astropy.table.Table` object.

        The user has to specify the name of the time column ``t_column`` to consider,
        the time ``format`` of this column, and its ``scale``.

        Parameters
        ----------
        table : `~astropy.table.Table`
            the astropy table from which the TimeMoc will be created
        t_column : str
            the name of the column in ``table`` containing the observation times
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)
        format : str, optional
            ``decimalyear`` by default. See `~astropy.time.Time` to know more about time scales and formats.
            ``format`` parameter is passed to the `~astropy.time.Time` time instantiations.
        scale : str, optional
            ``tdb`` by default. See `~astropy.time.Time` to know more about time scales and formats.
            ``scale`` parameter is passed to the `~astropy.time.Time` time instantiations.

        Returns
        -------
        moc : `~mocpy.tmoc.TimeMoc`
            the TimeMoc instance referring to the observations contained in ``table``

        """

        if not isinstance(t_column, str):
            raise TypeError

        moc = TimeMoc()

        times = Time(table[t_column],
                     format=kwargs.get('format', 'decimalyear'),
                     scale=kwargs.get('scale', 'tdb'))

        for time_start in times:
            moc.add_time(time_start)

        # degrade the TimeMoc to the order computer from ``delta_t``
        moc = moc.degrade_to_order(TimeMoc.time_resolution_to_order(delta_t))

        return moc

    @classmethod
    def from_csv_file(cls, path, delta_t=DEFAULT_OBSERVATION_TIME, **kwargs):
        """
        Load a TimeMoc from a CSV/ascii file containing two columns (t_min, t_max)

        Parameters
        ----------
        path : str
            CSV/ASCII path file
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)
        format : str, optional
            ``decimalyear`` by default. See `~astropy.time.Time` to know more about time scales and formats.
            ``format`` parameter is passed to the `~astropy.time.Time` time instantiations.
        scale : str, optional
            ``tdb`` by default. See `~astropy.time.Time` to know more about time scales and formats.
            ``scale`` parameter is passed to the `~astropy.time.Time` time instantiations.

        Returns
        -------
        moc : `~mocpy.tmoc.TimeMoc`
            the TimeMoc instance referring to the observations contained in the CSV file

        """

        moc = TimeMoc()

        table = Table.read(path)
        try:
            t_min = Time(table['t_min'],
                         format=kwargs.get('format', 'decimalyear'),
                         scale=kwargs.get('scale', 'tdb'))
            t_max = Time(table['t_max'],
                         format=kwargs.get('format', 'decimalyear'),
                         scale=kwargs.get('scale', 'tdb'))
        except KeyError:
            raise KeyError("'t_min' or 't_max' are not part of the table columns : {0}".format(table.info()))

        add_all_time_intervals = np.vectorize(moc.add_time_interval)
        add_all_time_intervals(time_start=t_min,
                               time_end=t_max)

        moc = moc.degrade_to_order(TimeMoc.time_resolution_to_order(delta_t))
        return moc

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

    def _process_degradation(self, another_moc, order_op):
        max_order = max(self.max_order, another_moc.max_order)
        if order_op > max_order:
            message = 'Requested time resolution for the operation cannot be applied.\n' \
                      'The TimeMoc object resulting from the operation is of time resolution {0} sec.'.format(
                TimeMoc.order_to_time_resolution(max_order).sec)
            warnings.warn(message, UserWarning)

        self_degradation = self.degrade_to_order(order_op)
        another_moc_degradation = another_moc.degrade_to_order(order_op)

        return self_degradation, another_moc_degradation

    def intersection(self, another_moc, delta_t=DEFAULT_OBSERVATION_TIME):
        order_op = TimeMoc.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMoc, self_degraded).intersection(moc_degraded)

    def union(self, another_moc, delta_t=DEFAULT_OBSERVATION_TIME):
        order_op = TimeMoc.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMoc, self_degraded).union(moc_degraded)

    def difference(self, another_moc, delta_t=DEFAULT_OBSERVATION_TIME):
        order_op = TimeMoc.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMoc, self_degraded).difference(moc_degraded)

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

    def filter_table(self, table, column_name, keep_inside=True, delta_t=DEFAULT_OBSERVATION_TIME, **kwargs):
        """
        Filter an astropy.table.Table to keep only rows inside (or outside) the tmoc instance

        Parameters
        ----------
        table : `~astropy.table.Table`
            the observations to filter
        column_name : str
            the name of the time column to consider for filtering the table.
        keep_inside : bool, optional
            describes if we want to keep rows inside or outside the time moc. True by default.
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations (i.e. the less precise order which
            is able to discriminate two observations separated by ``delta_t``)
        format : str, optional
            ``decimalyear`` by default. See `~astropy.time.Time` to know more about time scales and formats.
            ``format`` parameter is passed to the `~astropy.time.Time` time instantiations.
        scale : str, optional
            ``tdb`` by default. See `~astropy.time.Time` to know more about time scales and formats.
            ``scale`` parameter is passed to the `~astropy.time.Time` time instantiations.

        Returns
        -------
        result : `~astropy.table.Table`
            the filtered table (contains observations that are inside or outside the tmoc
            depending on the ``keep_inside`` parameter)

        """

        # the requested order for filtering the astropy observations table is more precise than the order
        # of the TimeMoc object
        current_max_order = self.max_order
        new_max_order = TimeMoc.time_resolution_to_order(delta_t)
        if new_max_order > current_max_order:
            message = 'Requested time resolution filtering cannot be applied.\n' \
                      'Filtering is applied with a time resolution of {0} sec.'.format(
                TimeMoc.order_to_time_resolution(current_max_order).sec)
            warnings.warn(message, UserWarning)

        rough_tmoc = self.degrade_to_order(new_max_order)

        time_arr = Time(table[column_name],
                        format=kwargs.get('format', 'decimalyear'),
                        scale=kwargs.get('scale', 'tdb'))
        pix_arr = (time_arr.jd * TimeMoc.DAY_MICRO_SEC)
        pix_arr = pix_arr.astype(int)

        intervals_arr = np.array(rough_tmoc._interval_set.intervals)
        inf_arr = np.vstack([pix_arr[i] >= intervals_arr[:, 0] for i in range(pix_arr.shape[0])])
        sup_arr = np.vstack([pix_arr[i] <= intervals_arr[:, 1] for i in range(pix_arr.shape[0])])

        if keep_inside:
            res = inf_arr & sup_arr
            filtered_rows = np.any(res, axis=1)
        else:
            res = ~inf_arr | ~sup_arr
            filtered_rows = np.all(res, axis=1)

        result = table[filtered_rows]
        return result

    def add_fits_header(self, tbhdu):
        """
        Add TIMESYS to the header of the fits describing the tmoc

        :param tbhdu: fits python object
        """
        tbhdu.header['TIMESYS'] = ('JD', 'ref system JD BARYCENTRIC TT, 1 usec level 29')

    @staticmethod
    def order_to_time_resolution(order):
        delta_t = TimeDelta(4**(29 - order) / 1e6, format='sec', scale='tdb')
        return delta_t

    @staticmethod
    def time_resolution_to_order(delta_time):
        order = 29 - int(np.log2(delta_time.sec * 1e6) / 2)
        return order

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
