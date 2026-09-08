import warnings
from copy import deepcopy

import numpy as np
from astropy.time import Time, TimeDelta

from .. import mocpy
from ..abstract_moc import AbstractMOC

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


DAY_MICRO_SEC = 86400000000.0


def times_to_microseconds(times):
    """
    Convert a `astropy.time.Time` into an array of integer microseconds since JD=0.

    This keeps the microsecond resolution required for `~mocpy.TimeMOC`.

    Parameters
    ----------
    times : `astropy.time.Time`
        Astropy observation times

    Returns
    -------
    `np.array`
        Time in microseconds

    Examples
    --------
    >>> from astropy.time import Time
    >>> from mocpy.tmoc import times_to_microseconds
    >>> time = Time("2026-08-15")
    >>> times_to_microseconds(time)
    np.uint64(212653512093460827)
    """
    times_jd = np.asarray(times.jd, dtype=np.uint64)
    times_us = np.asarray(
        (times - Time(times_jd, format="jd", scale="tcb")).jd * DAY_MICRO_SEC,
        dtype=np.uint64,
    )

    return times_jd * np.uint64(DAY_MICRO_SEC) + times_us


def microseconds_to_times(times_microseconds):
    """
    Convert an array of integer microseconds since JD=0, to an array of `astropy.time.Time`.

    Parameters
    ----------
    times_microseconds : `np.array`

    Returns
    -------
    `astropy.time.Time`

    Examples
    --------
    >>> from mocpy.tmoc import microseconds_to_times
    >>> time = microseconds_to_times(2e17)
    >>> time.iso
    '1625-08-24 07:33:20.000'
    """
    jd1 = np.asarray(times_microseconds // DAY_MICRO_SEC, dtype=np.float64)
    jd2 = np.asarray(
        (times_microseconds - jd1 * DAY_MICRO_SEC) / DAY_MICRO_SEC,
        dtype=np.float64,
    )

    return Time(val=jd1, val2=jd2, format="jd", scale="tcb")


class TimeMOC(AbstractMOC):
    """Multi-order time coverage class. Experimental."""

    # Maximum order of TimeMOCs
    # (do not remove since it may be used externally).
    MAX_ORDER = np.uint8(61)
    # Number of microseconds in a day
    DAY_MICRO_SEC = 86400000000.0
    # Default observation time : 30 min
    DEFAULT_OBSERVATION_TIME = TimeDelta(30 * 60, format="sec", scale="tcb")

    def __init__(self, store_index):
        """Is a Time Coverage (T-MOC).

        Args:
            store_index: index of the S-MOC in the rust-side storage
        """
        self.store_index = store_index

    @property
    def max_order(self):
        """Depth/order of the T-MOC."""
        depth = mocpy.get_tmoc_depth(self.store_index)
        return np.uint8(depth)

    @classmethod
    def n_cells(cls, depth):
        """Get the number of cells for a given depth.

        Parameters
        ----------
        depth : int
            The depth. It is comprised between 0 and `~mocpy.tmoc.TimeMOC.MAX_ORDER`

        Returns
        -------
        int
            The number of cells at the given order

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> TimeMOC.n_cells(0)
        2
        """
        if depth < 0 or depth > cls.MAX_ORDER:
            raise ValueError(
                f"The depth should be comprised between 0 and {cls.MAX_ORDER}, but {depth}"
                " was provided.",
            )
        return mocpy.n_cells_tmoc(depth)

    def to_time_ranges(self):
        """Return the time ranges this TimeMOC contains.

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> tmoc = TimeMOC.from_string("21/193392 22/386779-386780")
        >>> ranges = tmoc.to_time_ranges()
        >>> ranges.iso
        array([['2026-01-01 05:05:39.787', '2026-01-13 22:30:51.415'],
               ['2026-02-02 00:38:38.856', '2026-02-14 18:03:50.484']],
              dtype='<U23')
        """
        return microseconds_to_times(mocpy.to_ranges(self.store_index))

    @property
    def to_depth61_ranges(self):
        """Return the list of ranges this TimeMOC contains, in microsec since JD=0.

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> tmoc = TimeMOC.from_string("61/380000-390000 400000")
        >>> tmoc.to_depth61_ranges
        array([[380000, 390001],
               [400000, 400001]], dtype=uint64)
        """
        return mocpy.to_ranges(self.store_index)

    def degrade_to_order(self, new_order):
        """
        Degrade the MOC instance to a new, less precise, MOC.

        The maximum depth (i.e. the depth of the smallest Time cells that can be found in the MOC) of the
        degraded MOC is set to ``new_order``.

        Parameters
        ----------
        new_order : int
            Maximum depth of the output degraded MOC.

        Returns
        -------
        `~mocpy.TimeMOC`
            The degraded MOC.

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> tmoc = TimeMOC.from_string("20/0 100 500")
        >>> tmoc.max_order
        np.uint8(20)
        >>> tmoc.degrade_to_order(15)
        15/0 3 15
        """
        if new_order >= self.max_order:
            warnings.warn(
                "The new order is more precise than the current order, nothing done.",
                stacklevel=2,
            )
        index = mocpy.degrade(self.store_index, new_order)
        return TimeMOC(index)

    def refine_to_order(self, new_order):
        """Refine the order of the T-MOC instance to a more precise order.

        This is an in-place operation.

        Parameters
        ----------
        new_order : int
            New maximum order for this MOC.

        Returns
        -------
        `mocpy.TimeMOC`
            Returns itself, after in-place modification.

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> tmoc = TimeMOC.from_str("2/0")
        >>> tmoc
        2/0
        >>> tmoc.refine_to_order(3)
        2/0
        3/
        """
        if new_order <= self.max_order:
            warnings.warn(
                "'new_order' is less precise than the current max order. Nothing done.",
                stacklevel=2,
            )
        mocpy.refine(self.store_index, new_order)
        return self

    def to_order(self, new_order):
        """Create a new T-MOC with the new order.

        This is a convenience method for a quick change of order.
        Using 'degrade_to_order' and 'refine_to_order' depending on the situation is
        more efficient and avoids copying the MOC when it is not needed.

        Parameters
        ----------
        new_order : int
            The new order for the T-MOC. Can be either more or less precise than the
            current max_order of the T-MOC

        Returns
        -------
        `~mocpy.TimeMOC`
            A new T-MOC instance with the given max order.

        Examples
        --------
        >>> from mocpy import TimeMOC as TMOC
        >>> tmoc = TMOC.from_string("15/0-100")
        >>> tmoc.to_order(20)
        9/0
        10/2
        13/24
        15/100
        20/

        See Also
        --------
        degrade_to_order : to create a new less precise MOC
        refine_to_order : to change the order to a more precise one in place (no copy)
        """
        if new_order > self.max_order:
            moc_copy = deepcopy(self)
            return moc_copy.refine_to_order(new_order)
        if new_order < self.max_order:
            return self.degrade_to_order(new_order)
        return deepcopy(self)

    @classmethod
    def new_empty(cls, max_depth):
        """
        Create a new empty TimeMOC of given depth.

        Parameters
        ----------
        max_depth : int
            The resolution of the TimeMOC

        Returns
        -------
        `~mocpy.TimeMOC`

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> tmoc = TimeMOC.new_empty(10)
        >>> tmoc
        10/
        >>> tmoc.empty()
        True
        """
        index = mocpy.new_empty_tmoc(np.uint8(max_depth))
        return cls(index)

    @classmethod
    def from_depth61_ranges(cls, max_depth, ranges):
        """
        Create a TimeMOC from a set of Time ranges at order 61 (i.e. ranges of microseconds since JD=0).

        Parameters
        ----------
        max_depth : int
            The resolution of the TimeMOC
        ranges: `~numpy.ndarray`
                 a N x 2 numpy array representing the set of depth 61 ranges.

        Returns
        -------
        `~mocpy.TimeMOC`

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> import numpy as np
        >>> ranges = np.array([[0, 1000], [2000, 3000]], dtype=np.uint64)
        >>> tmoc = TimeMOC.from_depth61_ranges(61, ranges)
        >>> tmoc
        52/0 4
        53/2 10
        54/6 22
        55/14
        56/30 63 92
        57/125 186
        58/124 374
        61/
        """
        ranges = np.zeros((0, 2), dtype=np.uint64) if ranges is None else ranges

        if ranges.shape[1] != 2:
            raise ValueError(
                f"Expected a N x 2 numpy ndarray but second dimension is {ranges.shape[1]}",
            )

        if ranges.dtype is not np.uint64:
            ranges = ranges.astype(np.uint64)

        index = mocpy.from_time_ranges_array2(np.uint8(max_depth), ranges)
        return cls(index)

    @classmethod
    def from_times(cls, times, *, delta_t=DEFAULT_OBSERVATION_TIME, order=None):
        """
        Create a TimeMOC from a `astropy.time.Time`.

        Parameters
        ----------
        times : `astropy.time.Time`
            Astropy observation times
        delta_t : `astropy.time.TimeDelta`, optional
            The duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).
        order : int
            The order to use for the TimeMOC. If set, the `delta_t` will be ignored.

        Returns
        -------
        `~mocpy.TimeMOC`

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time, TimeDelta
        >>> delta = TimeDelta(1, format="jd")
        >>> times = Time(["2026-01-01", "2026-01-02", "2026-01-03"])
        >>> tmoc = TimeMOC.from_times(times, delta_t=delta, order=40)
        >>> tmoc.min_time.iso
        '2026-01-01 00:01:32.467'
        """
        times = times_to_microseconds(times)
        times = np.atleast_1d(times)

        if not order:
            order = TimeMOC.time_resolution_to_order(delta_t)
        store_index = mocpy.from_time_in_microsec_since_jd_origin(order, times)
        return cls(store_index)

    @classmethod
    def from_time_ranges(
        cls, min_times, max_times, *, delta_t=DEFAULT_OBSERVATION_TIME, order=None
    ):
        """
        Create a TimeMOC from a range defined by two `astropy.time.Time`.

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
        order : int
            The order to use for the TimeMOC. If set, the `delta_t` will be ignored.

        Returns
        -------
        `~mocpy.TimeMOC`

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min = Time(["2026-01-01", "2026-02-01"])
        >>> time_max = Time(["2026-01-20", "2026-02-20"])
        >>> tmoc = TimeMOC.from_time_ranges(time_min, time_max)
        >>> tmoc.contains(Time(["2026-01-15", "2027-01-01"]))
        array([ True, False])
        """
        if not order:
            # degrade the TimeMOC to the order computed from ``delta_t``
            order = TimeMOC.time_resolution_to_order(delta_t)

        min_times = times_to_microseconds(min_times)
        min_times = np.atleast_1d(min_times)

        max_times = times_to_microseconds(max_times)
        max_times = np.atleast_1d(max_times)

        if min_times.shape != max_times.shape:
            raise ValueError(
                f"Mismatch between min_times and max_times of shapes {min_times.shape} and {max_times.shape}",
            )

        store_index = mocpy.from_time_ranges_in_microsec_since_jd_origin(
            order,
            min_times,
            max_times,
        )
        return cls(store_index)

    @classmethod
    def from_time_ranges_approx(
        cls, min_times, max_times, *, delta_t=DEFAULT_OBSERVATION_TIME, order=None
    ):
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
        order : int
            The order to use for the TimeMOC. If set, the `delta_t` will be ignored.

        Returns
        -------
        `~mocpy.TimeMOC`

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min = Time("2026-01-01")
        >>> time_max = Time("2026-01-20")
        >>> tmoc = TimeMOC.from_time_ranges_approx(time_min, time_max)
        """
        if not order:
            # degrade the TimeMOC to the order computed from ``delta_t``
            order = TimeMOC.time_resolution_to_order(delta_t)

        min_times = np.asarray(min_times.jd)
        min_times = np.atleast_1d(min_times)

        max_times = np.asarray(max_times.jd)
        max_times = np.atleast_1d(max_times)
        if min_times.shape != max_times.shape:
            raise ValueError(
                f"Mismatch between min_times and max_times of shapes {min_times.shape} and {max_times.shape}",
            )

        store_index = mocpy.from_time_ranges(order, min_times, max_times)
        return cls(store_index)

    @classmethod
    def from_stmoc_space_fold(cls, smoc, stmoc):
        """
        Build a new T-MOC from the fold operation of the given ST-MOC by the given S-MOC.

        Parameters
        ----------
        smoc : `~mocpy.MOC`
            The Space-MOC to fold the ST-MOC with.
        stmoc : `~mocpy.STMOC`
            The Space-Time MOC the should be folded.

        Returns
        -------
        `~mocpy.TimeMOC`

        Examples
        --------
        >>> from mocpy import TimeMOC, MOC, STMOC
        >>> import astropy.units as u
        >>> from astropy.time import Time
        >>> smoc = MOC.from_cone(0*u.deg, 0*u.deg, radius=10*u.deg, max_depth=10)
        >>> stmoc = STMOC.from_spatial_coverages(Time("2000-01-01"), Time("2020-01-01"),
        ...                                      smoc, time_depth=40)
        >>> tmoc = TimeMOC.from_stmoc_space_fold(smoc, stmoc)
        >>> tmoc.min_time.iso
        '2000-01-01 00:01:14.142'
        """
        store_index = mocpy.project_on_stmoc_time_dim(
            smoc.store_index, stmoc.store_index
        )
        return cls(store_index)

    def _process_degradation(self, another_moc, order_op):
        """
        Degrade (down-sampling) self and ``another_moc`` to ``order_op`` order.

        Parameters
        ----------
        another_moc : `~mocpy.TimeMOC`
        order_op : int
            the order in which self and ``another_moc`` will be down-sampled to.

        Returns
        -------
        (`~mocpy.TimeMOC`, `~mocpy.TimeMOC`)
            self and ``another_moc`` degraded TimeMOCs

        """
        max_order = max(self.max_order, another_moc.max_order)
        if order_op > max_order:
            message = (
                "Requested time resolution for the operation cannot be applied.\n"
                f"The TimeMOC object resulting from the operation is of time resolution {TimeMOC.order_to_time_resolution(max_order).sec} sec."
            )
            warnings.warn(message, UserWarning, stacklevel=2)

        self_degradation = self.degrade_to_order(order_op)
        another_moc_degradation = another_moc.degrade_to_order(order_op)
        return self_degradation, another_moc_degradation

    def intersection_with_timeresolution(
        self, another_moc, *, delta_t=DEFAULT_OBSERVATION_TIME, order=None
    ):
        """
        Intersection between self and moc.

        ``delta_t`` gives the possibility to the user
        to set a time resolution for performing the tmoc intersection

        Parameters
        ----------
        another_moc : `~mocpy.TimeMOC`
            the TimeMOC used for performing the intersection with self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)
        order : int
            The order to use for the TimeMOC. If set, the `delta_t` will be ignored.

        Returns
        -------
        `~mocpy.TimeMOC`
            MOC object whose interval set corresponds to : self & ``moc``

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min1 = Time("2026-01-01")
        >>> time_max1 = Time("2026-01-10")
        >>> time_min2 = Time("2026-01-05")
        >>> time_max2 = Time("2026-01-15")
        >>> tmoc1 = TimeMOC.from_time_ranges(time_min1, time_max1)
        >>> tmoc2 = TimeMOC.from_time_ranges(time_min2, time_max2)
        >>> result = tmoc1.intersection_with_timeresolution(tmoc2)
        >>> result.to_time_ranges().iso
        array([['2026-01-04 23:45:57.301', '2026-01-10 00:15:48.998']],
              dtype='<U23')
        """
        if not order:
            order = TimeMOC.time_resolution_to_order(delta_t)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning, message="The new order*"
            )
            self_degraded, moc_degraded = self._process_degradation(another_moc, order)
        return super(TimeMOC, self_degraded).intersection(moc_degraded)

    def union_with_timeresolution(
        self, another_moc, *, delta_t=DEFAULT_OBSERVATION_TIME, order=None
    ):
        """
        Union between self and moc.

        ``delta_t`` gives the possibility to the user
        to set a time resolution for performing the tmoc union

        Parameters
        ----------
        another_moc : `~mocpy.TimeMOC`
            the TimeMOC to bind to self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)
        order : int
            The order to use for the TimeMOC. If set, the `delta_t` will be ignored.

        Returns
        -------
        `~mocpy.TimeMOC`
            MOC object whose interval set corresponds to : self | ``moc``

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min1 = Time("2026-01-01")
        >>> time_max1 = Time("2026-01-10")
        >>> time_min2 = Time("2026-01-15")
        >>> time_max2 = Time("2026-01-20")
        >>> tmoc1 = TimeMOC.from_time_ranges(time_min1, time_max1)
        >>> tmoc2 = TimeMOC.from_time_ranges(time_min2, time_max2)
        >>> result = tmoc1.union_with_timeresolution(tmoc2)
        >>> result.to_time_ranges().iso
        array([['2026-01-01 00:01:26.176', '2026-01-10 00:15:48.998'],
               ['2026-01-14 23:51:59.470', '2026-01-20 00:03:57.425']],
              dtype='<U23')
        """
        if not order:
            order = TimeMOC.time_resolution_to_order(delta_t)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning, message="The new order*"
            )
            self_degraded, moc_degraded = self._process_degradation(another_moc, order)
        return super(TimeMOC, self_degraded).union(moc_degraded)

    def difference_with_timeresolution(
        self, another_moc, *, delta_t=DEFAULT_OBSERVATION_TIME, order=None
    ):
        """Difference between self and another_moc.

        ``delta_t`` allows to set a time resolution to calculate the TimeMOC diff.

        Parameters
        ----------
        another_moc : `~mocpy.TimeMOC`
            the TimeMOC to substract from self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)
        order : int
            The order to use for the TimeMOC. If set, the `delta_t` will be ignored.

        Returns
        -------
        `~mocpy.TimeMOC`
            MOC object whose interval set corresponds to : self - ``moc``

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min1 = Time("2026-01-01")
        >>> time_max1 = Time("2026-01-20")
        >>> time_min2 = Time("2026-01-10")
        >>> time_max2 = Time("2026-01-15")
        >>> tmoc1 = TimeMOC.from_time_ranges(time_min1, time_max1)
        >>> tmoc2 = TimeMOC.from_time_ranges(time_min2, time_max2)
        >>> result = tmoc1.difference_with_timeresolution(tmoc2)
        >>> result.to_time_ranges().iso
        array([['2026-01-01 00:01:26.176', '2026-01-09 23:57:55.256'],
               ['2026-01-15 00:09:53.211', '2026-01-20 00:03:57.425']],
              dtype='<U23')
        """
        if not order:
            order = TimeMOC.time_resolution_to_order(delta_t)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning, message="The new order*"
            )
            self_degraded, moc_degraded = self._process_degradation(another_moc, order)
        return super(TimeMOC, self_degraded).difference(moc_degraded)

    @property
    def total_duration(self):
        """
        Get the total duration covered by the temporal moc.

        Returns
        -------
        `~astropy.time.TimeDelta`
            total duration of all the observation times of the tmoc

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min = Time("2026-01-01")
        >>> time_max = Time("2026-01-20")
        >>> tmoc = TimeMOC.from_time_ranges(time_min, time_max)
        >>> tmoc.total_duration.jd
        np.float64(19.001750565925924)
        """
        return TimeDelta(
            mocpy.ranges_sum(self.store_index) / 1e6,
            format="sec",
            scale="tcb",
        )

    @property
    def consistency(self):
        """
        Get a percentage of fill between the min and max time the moc is defined.

        A value near 0 shows a sparse temporal moc (i.e. the moc does not cover a lot
        of time and covers very distant times. A value near 1 means that the moc covers
        a lot of time without big pauses.

        Returns
        -------
        float
            fill percentage (between 0 and 1.)

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> times_min = Time(["2026-01-01", "2026-11-01"])
        >>> times_max = Time(["2026-01-20", "2026-11-30"])
        >>> tmoc = TimeMOC.from_time_ranges(times_min, times_max)
        >>> tmoc.consistency
        np.float64(0.14420062695924762)
        """
        return self.total_duration.jd / (self.max_time - self.min_time).jd

    @property
    def min_time(self):
        """
        Get the `~astropy.time.Time` time of the tmoc first observation.

        Returns
        -------
        `astropy.time.Time`
            time of the first observation

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min = Time("2026-01-01")
        >>> time_max = Time("2026-01-20")
        >>> tmoc = TimeMOC.from_time_ranges(time_min, time_max)
        >>> tmoc.min_time.iso
        '2026-01-01 00:01:26.176'
        """
        return microseconds_to_times(np.atleast_1d(self.min_index))[0]

    @property
    def max_time(self):
        """
        Get the `~astropy.time.Time` time of the tmoc last observation.

        Returns
        -------
        `~astropy.time.Time`
            time of the last observation

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min = Time("2026-01-01")
        >>> time_max = Time("2026-01-20")
        >>> tmoc = TimeMOC.from_time_ranges(time_min, time_max)
        >>> tmoc.max_time.iso
        '2026-01-20 00:03:57.425'
        """
        return microseconds_to_times(np.atleast_1d(self.max_index))[0]

    def contains(self, times, keep_inside=True):
        """
        Get a mask array (e.g. a numpy boolean array) of times being inside (or outside) the TMOC instance.

        Parameters
        ----------
        times : `astropy.time.Time`
            astropy times to check whether they are contained in the TMOC or not.
        keep_inside : bool, optional
            True by default. If so the filtered table contains only observations that are located the MOC.
            If ``keep_inside`` is False, the filtered table contains all observations lying outside the MOC.

        Returns
        -------
        `~numpy.array`
            A mask boolean array

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min = Time("2026-01-01")
        >>> time_max = Time("2026-01-20")
        >>> tmoc = TimeMOC.from_time_ranges(time_min, time_max)
        >>> tmoc.contains(Time(["2026-01-10", "2026-02-10"]))
        array([ True, False])
        """
        # the requested order for filtering the astropy observations table is more precise than the order
        # of the TimeMOC object
        pix_arr = times_to_microseconds(times)

        mask = mocpy.filter_time(self.store_index, pix_arr)

        if keep_inside:
            return mask
        return ~mask

    def contains_with_timeresolution(
        self, times, keep_inside=True, *, delta_t=DEFAULT_OBSERVATION_TIME, order=None
    ):
        """
        Get a mask array (e.g. a numpy boolean array) of times being inside (or outside) the TMOC instance.

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
        order : int
            The order to use for the TimeMOC. If set, the `delta_t` will be ignored.

        Returns
        -------
        `~numpy.array`
            A mask boolean array

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time, TimeDelta
        >>> time_min = Time("2026-01-01")
        >>> time_max = Time("2026-01-20")
        >>> delta = TimeDelta(2, format="jd")
        >>> tmoc = TimeMOC.from_time_ranges(time_min, time_max)
        >>> tmoc.contains_with_timeresolution(Time(["2026-01-20"]), delta_t=delta)
        array([ True])
        """
        # the requested order for filtering the astropy observations table is more precise than the order
        # of the TimeMOC object
        current_max_order = self.max_order
        new_max_order = order if order else TimeMOC.time_resolution_to_order(delta_t)
        if new_max_order > current_max_order:
            message = (
                "Requested time resolution filtering cannot be applied.\n"
                f"Filtering is applied with a time resolution of {TimeMOC.order_to_time_resolution(current_max_order).sec} sec."
            )
            warnings.warn(message, UserWarning, stacklevel=2)

        rough_tmoc = self.degrade_to_order(new_max_order)
        return rough_tmoc.contains(times, keep_inside)

    @staticmethod
    def order_to_time_resolution(order):
        """
        Convert an TimeMOC order to its equivalent time.

        Parameters
        ----------
        order : int
            order to convert

        Returns
        -------
        `~astropy.time.TimeDelta`
            time equivalent to ``order``

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> TimeMOC.order_to_time_resolution(61)
        <TimeDelta object: scale='tcb' format='sec' value=1e-06>
        >>> TimeMOC.order_to_time_resolution(0)
        <TimeDelta object: scale='tcb' format='sec' value=2305843009213.694>
        """
        return TimeDelta(2 ** (61 - order) / 1e6, format="sec", scale="tcb")

    @staticmethod
    def time_resolution_to_order(delta_time):
        """
        Convert a time resolution to a TimeMOC order.

        Parameters
        ----------
        delta_time : `~astropy.time.TimeDelta`
            time to convert

        Returns
        -------
        int
            The less precise order which is able to discriminate two observations separated by ``delta_time``.

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import TimeDelta
        >>> TimeMOC.time_resolution_to_order(TimeDelta(1, format="sec"))
        np.uint8(42)
        >>> TimeMOC.time_resolution_to_order(TimeDelta(30 * 60, format="sec"))
        np.uint8(31)
        """
        order = 61 - int(np.log2(delta_time.sec * 1e6))
        return np.uint8(order)

    def plot(
        self,
        *,
        title="TimeMOC",
        view=(None, None),
        figsize=(9.5, 5),
        cmap="Greens",
        **kwargs,
    ):
        """
        Plot the TimeMOC in a time window.

        This method uses interactive matplotlib: hover the plot to see the time in iso.

        Parameters
        ----------
        title : str, optional
            The title of the plot. Set to 'TimeMOC' by default.
        view : (`~astropy.time.Time`, `~astropy.time.Time`), optional
            Define the view window in which the observations are plotted.
            Set to (None, None) by default (i.e. all the observation time window
            is rendered).
        figsize: tuple[float], optional
            A tuple of two floats that will define the size of the figure.
        colors: tuple[str], optional
            A tuple of two matplotlib colors. The first will be used in the TimeMOC's
            holes while the second will represent the coverage of the TimeMOC.
        **kwargs:
            Any keyword argument that `~matplotlib.pyplot.imshow` accepts.

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time
        >>> time_min = Time(["2026-01-01", "2026-02-01"])
        >>> time_max = Time(["2026-01-20", "2026-02-20"])
        >>> tmoc = TimeMOC.from_time_ranges(time_min, time_max)
        >>> tmoc.plot() # doctest: +SKIP
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as err:
            raise ImportError("matplotlib is required to plot a TimeMOC.") from err

        if self.empty():
            import warnings

            warnings.warn("This time moc is empty", UserWarning, stacklevel=2)
            return

        plot_order = 30
        plotted_moc = (
            self.degrade_to_order(plot_order) if self.max_order > plot_order else self
        )

        min_jd = plotted_moc.min_time.jd if not view[0] else view[0].jd
        max_jd = plotted_moc.max_time.jd if not view[1] else view[1].jd

        if max_jd < min_jd:
            raise ValueError(
                f"Invalid selection: max_jd = {max_jd} must be > to min_jd = {min_jd}",
            )

        fig1 = plt.figure(figsize=figsize)
        ax = fig1.add_subplot(111)

        ax.set_xlabel("iso")
        ax.get_yaxis().set_visible(b=False)

        size = 2000
        delta = (max_jd - min_jd) / size
        min_jd_time = min_jd

        ax.set_xticks([0, size])
        ax.set_xticklabels(
            Time([min_jd_time, max_jd], format="jd", scale="tcb").iso,
            rotation=70,
        )

        y = np.zeros(size)
        for s_time_us, e_time_us in plotted_moc.to_time_ranges():
            s_index = int((s_time_us.jd - min_jd_time) / delta)
            e_index = int((e_time_us.jd - min_jd_time) / delta)
            y[s_index : (e_index + 1)] = 1.0

        # hack in case of full time mocs.
        if np.all(y):
            y[0] = 0

        z = np.tile(y, (int(size // 10), 1))

        plt.title(title)

        plt.imshow(z, interpolation="bilinear", cmap=cmap, **kwargs)

        label = ax.text(
            0,
            0,
            "",
            va="bottom",
            ha="left",
            fontsize=9,
            backgroundcolor="w",
            visible=False,
        )

        def on_mouse_motion(event):
            if event.inaxes is None or event.xdata is None or event.ydata is None:
                label.set_visible(False)
                fig1.canvas.draw_idle()
                return

            time = Time(event.xdata * delta + min_jd_time, format="jd", scale="tcb")

            tx = f"{time.iso}"
            label.set_position((event.xdata - 0.5, event.ydata - 0.5))
            label.set_text(tx)
            label.set_visible(True)
            fig1.canvas.draw_idle()

        fig1.canvas.mpl_connect("motion_notify_event", on_mouse_motion)

        plt.show()

    @classmethod
    def load(cls, path, format="fits"):  # noqa: A002
        """
        Load the Time MOC from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str or pathlib.Path
            The path to the file to load the MOC from.
        format : str, optional
            The format from which the MOC is loaded.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".

        Returns
        -------
        `~mocpy.TimeMOC`
        """
        path = str(path)
        if format == "fits":
            index = mocpy.time_moc_from_fits_file(path)
            return cls(index)
        if format == "ascii":
            index = mocpy.time_moc_from_ascii_file(path)
            return cls(index)
        if format == "json":
            index = mocpy.time_moc_from_json_file(path)
            return cls(index)
        formats = ("fits", "ascii", "json")
        raise ValueError(f"format should be one of {formats}")

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load MOC from raw bytes of a FITS file."""
        index = mocpy.time_moc_from_fits_raw_bytes(raw_bytes)
        return cls(index)

    @classmethod
    def from_string(cls, value, format="ascii"):  # noqa: A002
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

        Returns
        -------
        `~mocpy.TimeMOC`

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> tmoc = TimeMOC.from_string("5/0-10 8/200-300 305")
        >>> tmoc
        2/0
        3/7-8
        4/4 13
        5/10 25 36
        6/74
        8/300 305
        """
        if format == "ascii":
            index = mocpy.time_moc_from_ascii_str(value)
            return cls(index)
        if format == "json":
            index = mocpy.time_moc_from_json_str(value)
            return cls(index)
        formats = ("ascii", "json")
        raise ValueError(f"format should be one of {formats}")
