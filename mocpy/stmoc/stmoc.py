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
    Time-Spatial Coverage map class.

    """

    def __init__(self):
        self._index = None
        self._fits_column_name = 'PIXELS'

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

        if times.shape != lon.shape and lon.shape != lat.shape:
            raise ValueError("Times and positions must have the same length.")

        if times.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        result = cls()
        result._index = core.from_time_lonlat(times, time_depth, lon, lat, spatial_depth)
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
        if times.shape[1] != 2:
            raise ValueError("Times ranges must be provided. The shape of times must be (_, 2)")

        times = np.asarray(times.jd * 86400000000, dtype=np.uint64)
        ranges = core.project_on_second_dim(times, self._index)
        return MOC(IntervalSet(ranges, make_consistent=False))

    @property
    def _fits_header_keywords(self):
        t_depth, s_depth = core.ranges2d_depth(self._index)
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
        return core.ranges2d_to_uniq(self._index)