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

    def project_on_second_dimension(self, times):
        """
        Project the STMOC to its second dimension with a constraint on its first dimension.

        This will perform the union of all the spatial coverages lying in a set of time ranges.
        
        Parameters
        ----------
        times : `astropy.time.Time`
            The times ranges to project the STMOC onto.
        
        Returns
        -------
        result : `~mocpy.moc.MOC`
            The projeted Spatial Coverage map resulting from the union of all the
            spatial coverages lying in the set of time ranges given.
        """
        if times.ndim != 2 or times.shape[1] != 2:
            raise ValueError("Times ranges must be provided. The shape of times must be (_, 2)")

        times = np.asarray(times.jd * 86400000000, dtype=np.uint64)
        ranges = core.project_on_second_dim(times, self.__index)
        return MOC(IntervalSet(ranges, make_consistent=False))

    def project_on_first_dimension(self, moc):
        """
        Project the STMOC to its first dimension with a constraint
        on its second dimension.

        This will perform the union of all the time ranges whose associated
        spatial coverages lie in ``moc``.

        Parameters
        ----------
        moc : `mocpy.MOC`
            The spatial coverage to project the STMOC onto.

        Returns
        -------
        result : `~astropy.time.Time`
            A set of astropy defined time ranges.
        """
        # Time ranges in µsec
        time_ranges = core.project_on_first_dim(moc._intervals._intervals, self.__index)
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
