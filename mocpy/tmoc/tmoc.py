# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import warnings
import numpy as np

from astropy.time import Time, TimeDelta

from ..interval_set import IntervalSet
from ..abstract_moc import AbstractMOC

from .. import mocpy
from .. import utils

__author__ = "Matthieu Baumann, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "baumannmatthieu0@gmail.com, francois-xavier.pineau@astro.unistra.fr"

class TimeMOC(AbstractMOC):
    """Multi-order time coverage class. Experimental"""
    DAY_MICRO_SEC = 86400000000.
    # default observation time : 30 min
    DEFAULT_OBSERVATION_TIME = TimeDelta(30 * 60, format='sec', scale='tdb')

    # I introduced, but do not like, the double `make_consistent` (MOC + IntervalSet)
    # but `coverage_merge_time_intervals` is no more genric
    # and I can't remove `make_consistent` from `IntervalSet` without changing tests
    def __init__(self, interval_set=None, make_consistent=True, min_depth=None):
        """
        TimeMoc constructor.

        The merging step of the overlapping intervals is done here.

        Parameters
        ----------
        intervals : `~numpy.ndarray`
            a N x 2 numpy array representing the set of intervals.
        make_consistent : bool, optional
            True by default. Remove the overlapping intervals that makes
            a valid MOC (i.e. can be plot, serialized, manipulated).
        """
        super(TimeMOC, self).__init__(interval_set)

        if make_consistent:
            if min_depth is None:
                min_depth = -1

            min_depth = np.int8(min_depth)
            self._merge_intervals(min_depth)

    def _merge_intervals(self, min_depth):
        if not self.empty():
            self._interval_set._intervals = mocpy.coverage_merge_time_intervals(self._interval_set._intervals, min_depth)


    @property
    def max_order(self):
        """
        Depth of the smallest Time cells found in the MOC instance.
        """
        depth = mocpy.time_coverage_depth(self._interval_set._intervals)
        depth = np.uint8(depth)
        return depth

    def refine_to_order(self, min_depth):
        intervals = mocpy.coverage_merge_time_intervals(self._interval_set._intervals, min_depth)
        interval_set = IntervalSet(intervals, make_consistent=False)
        return TimeMOC(interval_set, make_consistent=False)

    def complement(self):
        """
        Returns the complement of the TimeMOC instance.

        Returns
        -------
        result : `~mocpy.moc.TimeMOC`
            The resulting TimeMOC.
        """
        intervals = mocpy.time_coverage_complement(self._interval_set._intervals)
        interval_set = IntervalSet(intervals, make_consistent=False)
        return TimeMOC(interval_set, make_consistent=False)

    def degrade_to_order(self, new_order):
        """
        Degrades the MOC instance to a new, less precise, MOC.

        The maximum depth (i.e. the depth of the smallest Time cells that can be found in the MOC) of the
        degraded MOC is set to ``new_order``.

        Parameters
        ----------
        new_order : int
            Maximum depth of the output degraded MOC.

        Returns
        -------
        moc : `~mocpy.tmoc.TimeMOC`
            The degraded MOC.
        """
        assert self._interval_set._intervals.shape[1] == 2
        intervals = mocpy.time_coverage_degrade(self._interval_set._intervals, new_order)
        return TimeMOC(IntervalSet(intervals, make_consistent=False), make_consistent=False)

    @classmethod
    def from_times(cls, times, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Create a TimeMOC from a `astropy.time.Time`

        Parameters
        ----------
        times : `astropy.time.Time`
            Astropy observation times
        delta_t : `astropy.time.TimeDelta`, optional
            The duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        time_moc : `~mocpy.tmoc.TimeMOC`
        """
        times = utils.times_to_microseconds(times)
        times = np.atleast_1d(times)

        times = times.reshape((times.shape[0], 1))
        intervals = np.hstack((times, times + np.uint64(1)))
        assert intervals.shape[1] == 2

        # degrade the TimeMoc to the order computed from ``delta_t``
        depth = TimeMOC.time_resolution_to_order(delta_t)
        tmoc = TimeMOC(IntervalSet(intervals))
        return tmoc.degrade_to_order(depth)

    @classmethod
    def from_time_ranges(cls, min_times, max_times, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Create a TimeMOC from a range defined by two `astropy.time.Time`

        Parameters
        ----------
        min_times : `astropy.time.Time`
            astropy times defining the left part of the intervals
        max_times : `astropy.time.Time`
            astropy times defining the right part of the intervals
        delta_t : `astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        time_moc : `~mocpy.tmoc.TimeMOC`
        """
        # degrade the TimeMoc to the order computed from ``delta_t``
        depth = TimeMOC.time_resolution_to_order(delta_t)
        
        min_times = utils.times_to_microseconds(min_times)
        min_times = np.atleast_1d(min_times)

        max_times = utils.times_to_microseconds(max_times)
        max_times = np.atleast_1d(max_times)

        assert min_times.shape == max_times.shape

        intervals = mocpy.from_time_ranges_in_microsec_since_jd_origin(min_times, max_times)

        tmoc = TimeMOC(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        return tmoc.degrade_to_order(depth)


    @classmethod
    def from_time_ranges_approx(cls, min_times, max_times, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Create a TimeMOC from a range defined by two `astropy.time.Time`.
        Uses the following approximation: simple take the JD time and multiply by the number of microseconds in a day.

        Parameters
        ----------
        min_times : `astropy.time.Time`
            astropy times defining the left part of the intervals
        max_times : `astropy.time.Time`
            astropy times defining the right part of the intervals
        delta_t : `astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        time_moc : `~mocpy.tmoc.TimeMOC`
        """
        # degrade the TimeMoc to the order computed from ``delta_t``
        depth = TimeMOC.time_resolution_to_order(delta_t)

        min_times = np.asarray(min_times.jd)
        min_times = np.atleast_1d(min_times)

        max_times = np.asarray(max_times.jd)
        max_times = np.atleast_1d(max_times)

        assert min_times.shape == max_times.shape

        intervals = mocpy.from_time_ranges(min_times, max_times)

        tmoc = TimeMOC(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        return tmoc.degrade_to_order(depth)

    def add_neighbours(self):
        """
        Add all the pixels at max order in the neighbourhood of the moc

        Returns
        -------
        tmoc : `~mocpy.tmoc.TimeMOC`
            self extended by one degree of neighbors.
        """
        time_delta = np.uint64(1) << (IntervalSet.TIME_MAX_ORDER - self.max_order)

        intervals = self._interval_set._intervals
        # WARN: astype gives the ownership/writeable of the array to python
        # It is necessary for writing the array
        # This will be removed as soon as this code is ported in rust
        intervals = intervals.astype(np.uint64)

        intervals[:, 0] = np.maximum(intervals[:, 0] - time_delta, np.uint64(0))
        intervals[:, 1] = np.minimum(intervals[:, 1] + time_delta, np.uint64((1 << 62) - 1))

        self._interval_set = IntervalSet(intervals)
        return self

    def remove_neighbours(self):
        """
        Remove all the pixels at max order located at the bound of the moc

        Returns
        -------
        tmoc : `~mocpy.tmoc.TimeMOC`
            self shrinked by one degree of neighbors.
        """
        time_delta = np.uint64(1) << (IntervalSet.TIME_MAX_ORDER - self.max_order)

        intervals = self._interval_set._intervals
        # WARN: astype gives the ownership/writeable of the array to python
        # It is necessary for writing the array
        # This will be removed as soon as this code is ported in rust
        intervals = intervals.astype(np.uint64)

        intervals[:, 0] = np.minimum(intervals[:, 0] + time_delta, np.uint64((1 << 62) - 1))
        intervals[:, 1] = np.maximum(intervals[:, 1] - time_delta, np.uint64(0))

        good_intervals = intervals[:, 1] > intervals[:, 0]

        self._interval_set = IntervalSet(intervals[good_intervals])
        return self

    def _process_degradation(self, another_moc, order_op):
        """
        Degrade (down-sampling) self and ``another_moc`` to ``order_op`` order

        Parameters
        ----------
        another_moc : `~mocpy.tmoc.TimeMoc`
        order_op : int
            the order in which self and ``another_moc`` will be down-sampled to.

        Returns
        -------
        result : (`~mocpy.tmoc.TimeMoc`, `~mocpy.tmoc.TimeMoc`)
            self and ``another_moc`` degraded TimeMocs

        """
        max_order = max(self.max_order, another_moc.max_order)
        if order_op > max_order:
            message = 'Requested time resolution for the operation cannot be applied.\n' \
                      'The TimeMoc object resulting from the operation is of time resolution {0} sec.'.format(
                TimeMOC.order_to_time_resolution(max_order).sec)
            warnings.warn(message, UserWarning)

        self_degradation = self.degrade_to_order(order_op)
        another_moc_degradation = another_moc.degrade_to_order(order_op)

        result = self_degradation, another_moc_degradation
        return result

    def intersection_with_timeresolution(self, another_moc, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Intersection between self and moc. ``delta_t`` gives the possibility to the user
        to set a time resolution for performing the tmoc intersection

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMOC used for performing the intersection with self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            MOC object whose interval set corresponds to : self & ``moc``

        """

        order_op = TimeMOC.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMOC, self_degraded).intersection(moc_degraded)

    def union_with_timeresolution(self, another_moc, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Union between self and moc. ``delta_t`` gives the possibility to the user
        to set a time resolution for performing the tmoc union

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMoc to bind to self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self | ``moc``

        """

        order_op = TimeMOC.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMOC, self_degraded).union(moc_degraded)

    def difference_with_timeresolution(self, another_moc, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Difference between self and moc. ``delta_t`` gives the possibility to the user
        to set a time resolution for performing the tmoc diff

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMoc to substract from self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self - ``moc``

        """

        order_op = TimeMOC.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMOC, self_degraded).difference(moc_degraded)

    @property
    def _fits_header_keywords(self):
        return {
            'PIXTYPE': 'HEALPIX',
            'ORDERING': 'NUNIQ',
            'TIMESYS': ('JD', 'ref system JD BARYCENTRIC TT, 1 usec level 29'),
            'MOCORDER': self.max_order,
            'MOCTOOL': 'MOCPy'
        }
    
    @property
    def _fits_format(self):
        depth = self.max_order
        if depth <= 13:
            fits_format = '1J'
        else:
            fits_format = '1K'
        return fits_format

    @property
    def total_duration(self):
        """
        Get the total duration covered by the temporal moc

        Returns
        -------
        duration : `~astropy.time.TimeDelta`
            total duration of all the observation times of the tmoc
            total duration of all the observation times of the tmoc

        """

        if self._interval_set.empty():
            return 0

        total_time_us = 0
        # The interval set is checked for consistency before looping over all the intervals
        for (start_time, stop_time) in self._interval_set._intervals:
            total_time_us = total_time_us + (stop_time - start_time)

        duration = TimeDelta(total_time_us / 1e6, format='sec', scale='tdb')
        return duration

    @property
    def consistency(self):
        """
        Get a percentage of fill between the min and max time the moc is defined.

        A value near 0 shows a sparse temporal moc (i.e. the moc does not cover a lot
        of time and covers very distant times. A value near 1 means that the moc covers
        a lot of time without big pauses.

        Returns
        -------
        result : float
            fill percentage (between 0 and 1.)

        """

        result = self.total_duration.jd / (self.max_time - self.min_time).jd
        return result

    @property
    def min_time(self):
        """
        Get the `~astropy.time.Time` time of the tmoc first observation

        Returns
        -------
        min_time : `astropy.time.Time`
            time of the first observation

        """

        # min_time = Time(self._interval_set.min / TimeMOC.DAY_MICRO_SEC, format='jd', scale='tdb')
        min_time = utils.microseconds_to_times(self._interval_set.min)
        return min_time

    @property
    def max_time(self):
        """
        Get the `~astropy.time.Time` time of the tmoc last observation

        Returns
        -------
        max_time : `~astropy.time.Time`
            time of the last observation

        """

        max_time = utils.microseconds_to_times(self._interval_set.max)
        return max_time

    def contains(self, times, keep_inside=True, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Get a mask array (e.g. a numpy boolean array) of times being inside (or outside) the
        TMOC instance.

        Parameters
        ----------
        times : `astropy.time.Time`
            astropy times to check whether they are contained in the TMOC or not.
        keep_inside : bool, optional
            True by default. If so the filtered table contains only observations that are located the MOC.
            If ``keep_inside`` is False, the filtered table contains all observations lying outside the MOC.
        delta_t : `astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        array : `~numpy.darray`
            A mask boolean array
        """
        # the requested order for filtering the astropy observations table is more precise than the order
        # of the TimeMoc object
        current_max_order = self.max_order
        new_max_order = TimeMOC.time_resolution_to_order(delta_t)
        if new_max_order > current_max_order:
            message = 'Requested time resolution filtering cannot be applied.\n' \
                      'Filtering is applied with a time resolution of {0} sec.'.format(
                TimeMOC.order_to_time_resolution(current_max_order).sec)
            warnings.warn(message, UserWarning)

        rough_tmoc = self.degrade_to_order(new_max_order)

        #pix_arr = (times.jd * TimeMOC.DAY_MICRO_SEC)
        #pix_arr = pix_arr.astype(np.uint64)
        pix_arr = utils.times_to_microseconds(times)

        intervals = rough_tmoc._interval_set._intervals
        inf_arr = np.vstack([pix_arr[i] >= intervals[:, 0] for i in range(pix_arr.shape[0])])
        sup_arr = np.vstack([pix_arr[i] <= intervals[:, 1] for i in range(pix_arr.shape[0])])

        if keep_inside:
            res = inf_arr & sup_arr
            filtered_rows = np.any(res, axis=1)
        else:
            res = ~inf_arr | ~sup_arr
            filtered_rows = np.all(res, axis=1)

        return filtered_rows

    @staticmethod
    def order_to_time_resolution(order):
        """
        Convert an TimeMoc order to its equivalent time

        Parameters
        ----------
        order : int
            order to convert

        Returns
        -------
        delta_t : `~astropy.time.TimeDelta`
            time equivalent to ``order``

        """

        delta_t = TimeDelta(2**(61 - order) / 1e6, format='sec', scale='tdb')
        return delta_t

    @staticmethod
    def time_resolution_to_order(delta_time):
        """
        Convert a time resolution to a TimeMoc order.

        Parameters
        ----------
        delta_time : `~astropy.time.TimeDelta`
            time to convert

        Returns
        -------
        order : int
            The less precise order which is able to discriminate two observations separated by ``delta_time``.

        """

        order = 61 - int(np.log2(delta_time.sec * 1e6))
        return np.uint8(order)

    def plot(self, title='TimeMoc', view=(None, None), figsize=(9.5, 5), **kwargs):
        """
        Plot the TimeMoc in a time window.

        This method uses interactive matplotlib. The user can move its mouse through the plot to see the
        time (at the mouse position).

        Parameters
        ----------
        title : str, optional
            The title of the plot. Set to 'TimeMoc' by default.
        view : (`~astropy.time.Time`, `~astropy.time.Time`), optional
            Define the view window in which the observations are plotted. Set to (None, None) by default (i.e.
            all the observation time window is rendered).

        """
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.pyplot as plt

        if self._interval_set.empty():
            import warnings
            warnings.warn('This time moc is empty', UserWarning)
            return

        plot_order = 30
        if self.max_order > plot_order:
            plotted_moc = self.degrade_to_order(plot_order)
        else:
            plotted_moc = self

        min_jd = plotted_moc.min_time.jd if not view[0] else view[0].jd
        max_jd = plotted_moc.max_time.jd if not view[1] else view[1].jd

        
        if max_jd < min_jd:
            raise ValueError("Invalid selection: max_jd = {0} must be > to min_jd = {1}".format(max_jd, min_jd))

        fig1 = plt.figure(figsize=figsize)
        ax = fig1.add_subplot(111)

        ax.set_xlabel('iso')
        ax.get_yaxis().set_visible(False)

        size = 2000
        delta = (max_jd - min_jd) / size
        min_jd_time = min_jd

        ax.set_xticks([0, size])
        ax.set_xticklabels(Time([min_jd_time, max_jd], format='jd', scale='tdb').iso, rotation=70)

        y = np.zeros(size)
        for (s_time_us, e_time_us) in plotted_moc._interval_set._intervals:
            s_index = int((utils.microseconds_to_times(s_time_us).jd - min_jd_time) / delta)
            e_index = int((utils.microseconds_to_times(e_time_us).jd - min_jd_time) / delta)
            y[s_index:(e_index+1)] = 1.0

        # hack in case of full time mocs.
        if np.all(y):
            y[0] = 0

        z = np.tile(y, (int(size//10), 1))

        plt.title(title)

        color_map = LinearSegmentedColormap.from_list('w2r', ['#fffff0', '#aa0000'])
        color_map.set_under('w')
        color_map.set_bad('gray')

        plt.imshow(z, interpolation='bilinear', **kwargs)

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

    def save(self, path, format='fits', overwrite=False):
        """
        Writes the Time MOC to a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str
            The path to the file to save the MOC in.
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        overwrite : bool, optional
            If the file already exists and you want to overwrite it, then set the  ``overwrite`` keyword.
            Default to False.
        """
        import os
        file_exists = os.path.isfile(path)
        
        if file_exists and not overwrite:
            raise OSError('File {} already exists! Set ``overwrite`` to '
                          'True if you want to replace it.'.format(path))        
        
        if format == 'fits':
            mocpy.time_moc_to_fits_file(self.max_order, self._interval_set._intervals, path)
        elif format == 'ascii':
            mocpy.time_moc_to_ascii_file(self.max_order, self._interval_set._intervals, path)
        elif format == 'json':
            mocpy.time_moc_to_json_file(self.max_order, self._interval_set._intervals, path)
        else:
            formats = ('fits', 'ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))

    @classmethod
    def load(cls, path, format='fits'):
        """
        Load the Time MOC from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str
            The path to the file to load the MOC from.
        format : str, optional
            The format from which the MOC is loaded.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        """
        if format == 'fits':
            intervals = mocpy.time_moc_from_fits_file(path)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        elif format == 'ascii':
            intervals = mocpy.time_moc_from_ascii_file(path)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        elif format == 'json':
            intervals = mocpy.time_moc_from_json_file(path)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        else:
            formats = ('fits', 'ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))

    def to_string(self, format='ascii'):
        """
        Writes the Time MOC into a string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == 'ascii':
            return mocpy.time_moc_to_ascii_str(self.max_order, self._interval_set._intervals)
        elif format == 'json':
            return mocpy.time_moc_to_json_str(self.max_order, self._interval_set._intervals)
        else:
            formats = ('ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))

    @classmethod
    def from_string(cls, value, format='ascii'):
        """
        Deserialize the Time MOC from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        WARNING: the serialization must be strict, i.e. **must not** contain overlapping elements

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == 'ascii':
            intervals = mocpy.time_moc_from_ascii_str(value)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        elif format == 'json':
            intervals = mocpy.time_moc_from_json_str(value)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        else:
            formats = ('ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))