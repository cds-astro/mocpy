# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
from urllib.parse import urlencode
from io import BytesIO

from astropy.utils.data import download_file
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import ICRS, Galactic, BaseCoordinateFrame
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.time import Time

import cdshealpix
from .. import core

from .. import MOC, serializer
from ..interval_set import IntervalSet

__author__ = "Thomas Boch, Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, matthieu.baumann@astro.unistra.fr"


class STMOC(serializer.IO):
    """
    Time-Spatial Coverage class.

    """

    def __init__(self):
        self.__index = core.create_2d_coverage()
        self._fits_column_name = 'PIXELS'

    def __del__(self):
        core.drop_2d_coverage(self.__index)

    def __eq__(self, other):
        return core.coverage_2d_equality_check(self.__index, other.__index)

    @property
    def max_depth(self):
        return core.coverage_2d_depth(self.__index)

    @property
    def max_time(self):
        return Time(core.coverage_2d_max_time(self.__index), format='jd', scale='tdb')
    
    @property
    def min_time(self):
        return Time(core.coverage_2d_min_time(self.__index), format='jd', scale='tdb')

    def is_empty(self):
        """
        Check whether the Space-Time coverage is empty
        """
        return core.coverage_2d_is_empty(self.__index)

    @classmethod
    def from_times_positions(cls, times, time_depth, lon, lat, spatial_depth):
        """
        Creates a 2D Coverage from a set of times and positions associated to each time.

        - Its first dimension refers to `astropy.time.Time` times.
        - Its second dimension refers to lon, lat `astropy.units.Quantity` positions.
        
        Parameters
        ----------
        time : `astropy.time.Time`
            The times of each sky coordinates.
        time_depth : int
            Time depth.
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed at a specific time.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed at a specific time.
        spatial_depth : int
            Spatial depth.
        
        Returns
        -------
        result : `~mocpy.stmoc.STMOC`
            The resulting Spatial-Time Coverage map.
        """
        times = times.jd.astype(np.float64)
        lon = lon.to_value('rad').astype(np.float64)
        lat = lat.to_value('rad').astype(np.float64)

        if times.shape != lon.shape or lon.shape != lat.shape:
            raise ValueError("Times and positions must have the same length.")

        if times.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        result = cls()
        core.from_time_lonlat(result.__index, times, time_depth, lon, lat, spatial_depth)
        return result

    @classmethod
    def from_time_ranges_positions(cls, times_start, times_end, lon, lat, time_depth=29, spatial_depth=29):
        """
        Creates a 2D Coverage from a set of times and positions associated to each time.

        - Its first dimension refers to `astropy.time.Time` times.
        - Its second dimension refers to lon, lat `astropy.units.Quantity` positions.

        Parameters
        ----------
        times_start : `astropy.time.Time`
            The starting times of each observations.
        times_end : `astropy.time.Time`
            The ending times of each observations.
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed at a specific time.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed at a specific time.
        time_depth : int, optional
            Time depth. By default, the time resolution chosen is 1µs.
        spatial_depth : int, optional
            Spatial depth. By default, the space resolution chosen is 393.2μas.

        Returns
        -------
        result : `~mocpy.stmoc.STMOC`
            The resulting Spatial-Time Coverage map.
        """
        times_start = times_start.jd.astype(np.float64)
        times_end = times_end.jd.astype(np.float64)

        lon = lon.to_value('rad').astype(np.float64)
        lat = lat.to_value('rad').astype(np.float64)

        if times_start.shape != lon.shape or lon.shape != lat.shape or times_start.shape != times_end.shape:
            raise ValueError("Times and positions must have the same length.")

        if times_start.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        result = cls()
        core.from_time_ranges_lonlat(result.__index, times_start, times_end, time_depth, lon, lat, spatial_depth)
        return result

    @classmethod
    def from_spatial_coverages(cls, times_start, times_end, spatial_coverages, time_depth=29):
        """
        Creates a 2D Coverage from a set of time ranges and spatial coverages associated to them.

        Parameters
        ----------
        times_start : `astropy.time.Time`
            The starting times of each observations.
        times_end : `astropy.time.Time`
            The ending times of each observations.
        spatial_coverages : list
            List of `mocpy.MOC` spatial coverage objects.
        time_depth : int, optional
            Time depth. By default, the time resolution chosen is 1µs.
        spatial_depth : int, optional
            Spatial depth. By default, the space resolution chosen is 393.2μas.

        Returns
        -------
        result : `~mocpy.stmoc.STMOC`
            The resulting Spatial-Time Coverage map.
        """
        times_start = times_start.jd.astype(np.float64)
        times_end = times_end.jd.astype(np.float64)

        if times_start.shape != times_end.shape or times_start.shape[0] != len(spatial_coverages):
            raise ValueError("Time ranges and spatial coverages must have the same length")

        if times_start.ndim != 1:
            raise ValueError("Times and spatial coverages must be 1D arrays")

        result = cls()
        spatial_coverages = [spatial_coverage._interval_set._intervals for spatial_coverage in spatial_coverages]

        core.from_time_ranges_spatial_coverages(result.__index, times_start, times_end, time_depth, spatial_coverages)
        return result

    def query_by_time(self, times):
        """
        Query the ST-MOC by time ranges.

        This will perform the union of all the spatial coverages lying in a set of time ranges.

        Parameters
        ----------
        times : `astropy.time.Time`
            Time ranges. Must be a Nx2 shaped astropy time array.
        
        Returns
        -------
        result : `~mocpy.moc.MOC`
            The spatial coverage being observed within the input time ranges
        """
        if times.ndim != 2 or times.shape[1] != 2:
            raise ValueError("Times ranges must be provided. The shape of times must be (_, 2)")

        times = np.asarray(times.jd * 86400000000, dtype=np.uint64)
        ranges = core.project_on_second_dim(times, self.__index)
        return MOC(IntervalSet(ranges, make_consistent=False))

    def query_by_space(self, spatial_coverage):
        """
        Query the ST-MOC by space coverage.

        This will perform the union of all the time ranges whose associated
        spatial coverages lie in ``moc``.

        Parameters
        ----------
        spatial_coverage : `mocpy.MOC`
            The spatial coverage.

        Returns
        -------
        result : `~astropy.time.Time`
            The time ranges observing in the ``spatial_coverage``
        """
        # Time ranges in µsec
        time_ranges = core.project_on_first_dim(spatial_coverage._interval_set._intervals, self.__index)
        return Time(time_ranges / 86400000000, format='jd', scale='tdb')

    def union(self, other):
        """
        Union of two Space-Time coverages.

        The union is first performed along the Time axis of the
        coverage. The union of space axis coverages is done on time
        ranges that overlap.

        Parameters
        ----------
        other : `mocpy.STMOC`
            The Space-Time coverage to perform the union with.

        Returns
        -------
        result : `mocpy.STMOC`
            A new Space-Time coverage being the union of `self`
            with `other`.
        """
        result = STMOC()
        core.coverage_2d_union(result.__index, self.__index, other.__index)
        return result

    def intersection(self, other):
        """
        Intersection of two Space-Time coverages.

        The intersection is first performed along the Time axis of the
        coverage. The intersection of space axis coverages is done on time
        ranges that overlap.

        Parameters
        ----------
        other : `mocpy.STMOC`
            The Space-Time coverage to perform the intersection with.

        Returns
        -------
        result : `mocpy.STMOC`
            A new Space-Time coverage being the intersection of `self`
            with `other`.
        """
        result = STMOC()
        core.coverage_2d_intersection(result.__index, self.__index, other.__index)
        return result

    def difference(self, other):
        """
        Difference of two Space-Time coverages.

        The difference is first performed along the Time axis of the
        coverage. The difference of space axis coverages is done on time
        ranges that overlap.

        Parameters
        ----------
        other : `mocpy.STMOC`
            The Space-Time coverage to perform the difference with.

        Returns
        -------
        result : `mocpy.STMOC`
            A new Space-Time coverage being the difference of `self`
            with `other`.
        """
        result = STMOC()
        core.coverage_2d_difference(result.__index, self.__index, other.__index)
        return result

    def contains(self, times, lon, lat, inside=True):
        """
        Returns a boolean mask array of the (times, positions) lying inside (or outside) the Space-Time coverage.

        Parameters
        ----------
        times : `astropy.time.Time`
            The times of each sky coordinates.
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed at a specific time.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed at a specific time.
        inside : bool, optional
            True by default. The returned mask array has true values for (time, position)
            lying inside the Space-Time coverage.

        Raises
        ------
        ValueError : If `times`, `lon` and `lat` do not have the same length.

        Returns
        -------
        array : `~np.ndarray`
            A mask boolean array
        """
        times = times.jd.astype(np.float64)
        lon = lon.to_value('rad').astype(np.float64)
        lat = lat.to_value('rad').astype(np.float64)

        if times.shape != lon.shape or lon.shape != lat.shape:
            raise ValueError("Times and positions must have the same length.")

        if times.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        result = core.coverage_2d_contains(self.__index, times, lon, lat)

        if not inside:
            result = ~result

        return result

    @property
    def _fits_header_keywords(self):
        t_depth, s_depth = core.coverage_2d_depth(self.__index)
        return {
            'MOC': 'TIME.SPACE',
            'ORDERING': 'RANGE29',
            # Max depth along the first dimension
            'MOCORDER': str(t_depth),
            # Max depth along the second dimension
            'MOCORD_1': str(s_depth),
            'COORDSYS': 'C',
            'TIMESYS': 'JD',
            'MOCTOOL': 'MOCPy'
        }

    @property
    def _fits_format(self):
        return '1K'

    def _uniq_format(self):
        return core.coverage_2d_to_fits(self.__index)

    @classmethod
    def deserialization(cls, hdulist):
        """
        Deserialization of an hdulist to a Time-Space coverage
        """
        # The binary HDU table contains all the data
        header = hdulist[1].header

        bin_HDU_table = hdulist[1]
        # Some FITS file may not have a name column to access
        # the data. So we rename it as the `PIXELS` columns.
        if len(bin_HDU_table.columns) != 1:
            raise AttributeError('The binary HDU table of {0} must contain only one'
                'column where the data are stored.')

        if bin_HDU_table.columns[0].name is None:
            bin_HDU_table.columns[0].name = 'PIXELS'

        key = bin_HDU_table.columns[0].name

        # Retrieve the spatial and time order
        # FITS header can be of two different format for expressing
        # these orders:
        # 1. TORDER key refers the max depth in the time dimension. MOCORDER then
        #    refers the spatial max depth.
        # 2. MOCORDER refers the max depth along the 1st dimension whereas
        #    MOCORD_1 refers to the max depth along the 2nd dimension.
        if header.get('TORDER') is None:
            # We are in the 2. case
            first_dim_depth = header.get('MOCORDER')
            second_dim_depth = header.get('MOCORD_1')
        else:
            # We are in the 1. case
            first_dim_depth = header.get('TORDER')
            second_dim_depth = header.get('MOCORDER')

        result = cls()
        core.coverage_2d_from_fits(
            result.__index,
            bin_HDU_table.data[key].astype(np.int64)
        )
        return result

    @classmethod
    def from_fits(cls, filename):
        """
        Loads a STMOC from a FITS file.

        Parameters
        ----------
        filename : str
            The path to the FITS file.

        Returns
        -------
        result : `~mocpy.moc.STMOC`
            The resulting STMOC.
        """
        # Open the FITS file
        hdulist = fits.open(filename)
        return STMOC.deserialization(hdulist)
