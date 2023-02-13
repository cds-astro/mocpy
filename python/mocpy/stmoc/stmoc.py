# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import numpy as np

from .. import MOC, mocpy
from ..tmoc import TimeMOC, microseconds_to_times, times_to_microseconds
from ..abstract_moc import AbstractMOC


__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


class STMOC(AbstractMOC):
    """Time-Spatial Coverage class."""

    __create_key = object()

    def __init__(self, create_key, store_index):
        """Is a Time-Spatial Coverage (ST-MOC).

        Args:
            create_key: Object ensure __init__ is called by super-class/class-methods only
            store_index: index of the ST-MOC in the rust-side storage
        """
        super(STMOC, self).__init__(
            AbstractMOC._create_key, STMOC.__create_key, store_index
        )
        assert (
            create_key == STMOC.__create_key
        ), "ST-MOC instantiation is only allowed by class or super-class methods"

    def __del__(self):
        """Erase STMOC."""
        super(STMOC, self).__del__()

    def __eq__(self, other):
        """Assert equality between MOCs."""
        return super(STMOC, self).__eq__(other)

    @property
    def max_depth(self):
        """Return max depth of MOC."""
        return mocpy.coverage_2d_depth(self._store_index)

    @property
    def max_time(self):
        """Return STMOC max time."""
        return microseconds_to_times(mocpy.coverage_2d_max_time(self._store_index))

    @property
    def min_time(self):
        """Return STMOC min time."""
        return microseconds_to_times(mocpy.coverage_2d_min_time(self._store_index))

    def is_empty(self):
        """Check whether the Space-Time coverage is empty."""
        return mocpy.is_empty(self._store_index)

    @classmethod
    def from_times_positions(cls, times, time_depth, lon, lat, spatial_depth):
        """
        Create a 2D Coverage from a set of times and positions associated to each time.

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
        times = times_to_microseconds(times)
        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)

        if times.shape != lon.shape or lon.shape != lat.shape:
            raise ValueError("Times and positions must have the same length.")

        if times.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        index = mocpy.from_time_lonlat(times, time_depth, lon, lat, spatial_depth)
        result = cls(cls.__create_key, index)

        return result

    @classmethod
    def from_time_ranges_positions(
        cls, times_start, times_end, lon, lat, time_depth=61, spatial_depth=29
    ):
        """
        Create a 2D Coverage from a set of times and positions associated to each time.

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
        # times_start = times_start.jd.astype(np.float64)
        # times_end = times_end.jd.astype(np.float64)

        times_start = times_to_microseconds(times_start)
        times_end = times_to_microseconds(times_end)

        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)

        if (
            times_start.shape != lon.shape
            or lon.shape != lat.shape
            or times_start.shape != times_end.shape
        ):
            raise ValueError("Times and positions must have the same length.")

        if times_start.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        index = mocpy.from_time_ranges_lonlat(
            times_start, times_end, time_depth, lon, lat, spatial_depth
        )
        result = cls(cls.__create_key, index)

        return result

    @classmethod
    def from_spatial_coverages(
        cls, times_start, times_end, spatial_coverages, time_depth=61
    ):
        """
        Create a 2D Coverage from a set of time ranges and spatial coverages associated to them.

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

        Returns
        -------
        result : `~mocpy.stmoc.STMOC`
            The resulting Spatial-Time Coverage map.
        """
        # times_start = times_start.jd.astype(np.float64)
        # times_end = times_end.jd.astype(np.float64)
        times_start = times_to_microseconds(times_start)
        times_end = times_to_microseconds(times_end)

        if times_start.shape != times_end.shape or times_start.shape[0] != len(
            spatial_coverages
        ):
            raise ValueError(
                "Time ranges and spatial coverages must have the same length"
            )

        if times_start.ndim != 1:
            raise ValueError("Times and spatial coverages must be 1D arrays")


        spatial_coverages_indices = np.fromiter(
            (arg._store_index for arg in spatial_coverages), dtype=AbstractMOC.store_index_dtype()
        )
        index = mocpy.from_time_ranges_spatial_coverages(
            times_start, times_end, time_depth, spatial_coverages_indices
        )
        result = cls(cls.__create_key, index)

        return result

    def query_by_time(self, tmoc):
        """
        Query the ST-MOC by time T-MOC.

        This will perform the union of all the spatial coverages lying in a set of time ranges.

        Parameters
        ----------
        tmoc : ~mocpy.tmoc.TimeMOC``
            Time ranges. Must be a Nx2 shaped astropy time array.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The spatial coverage being observed within the input time ranges
        """
        # if times.ndim != 2 or times.shape[1] != 2:
        #    raise ValueError(
        #        "Times ranges must be provided. The shape of times must be (_, 2)"
        #    )

        # times = np.asarray(times.jd * 86400000000, dtype=np.uint64)
        # times = utils.times_to_microseconds(times)
        # ranges = mocpy.project_on_second_dim(times, self._store_index)
        # return MOC(IntervalSet(ranges, make_consistent=False))
        # store_index = from_stmoc_time_fold(cls, tmoc, stmoc):
        return MOC.from_stmoc_time_fold(tmoc, self)

    def query_by_space(self, smoc):
        """
        Query the ST-MOC by space coverage.

        This will perform the union of all the time ranges whose associated
        spatial coverages lie in ``moc``.

        Parameters
        ----------
        smoc : `~mocpy.moc.MOC`
            The spatial coverage.

        Returns
        -------
        result : `~mocpy.tmoc.TimeMOC`
            The time ranges observing in the ``spatial_coverage``
        """
        return TimeMOC.from_stmoc_space_fold(smoc, self)

    def contains(self, times, lon, lat, inside=True):
        """
        Return a boolean mask array of the (times, positions) lying inside (or outside) the Space-Time coverage.

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
        # times = times.jd.astype(np.float64)
        times = times_to_microseconds(times)
        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)

        if times.shape != lon.shape or lon.shape != lat.shape:
            raise ValueError("Times and positions must have the same length.")

        if times.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        result = mocpy.coverage_2d_contains(self._store_index, times, lon, lat)

        if not inside:
            result = ~result

        return result

    @classmethod
    def load(cls, path, format="fits"):
        """
        Load the Spatial MOC from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str or pathlib.Path
            The path to the file to load the MOC from.
        format : str, optional
            The format from which the MOC is loaded.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        """
        path = str(path)
        if format == "fits":
            index = mocpy.coverage_2d_from_fits_file(path)
            return cls(cls.__create_key, index)
        elif format == "ascii":
            index = mocpy.coverage_2d_from_ascii_file(path)
            return cls(cls.__create_key, index)
        elif format == "json":
            index = mocpy.coverage_2d_from_json_file(path)
            return cls(cls.__create_key, index)
        else:
            formats = ("fits", "ascii", "json")
            raise ValueError("format should be one of %s" % (str(formats)))

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load MOC from raw bytes of a FITS file."""
        index = mocpy.stmoc_from_fits_raw_bytes(raw_bytes)
        return cls(cls.__create_key, index)

    @classmethod
    def from_string(cls, value, format="ascii"):
        """
        Deserialize the Spatial MOC from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == "ascii":
            index = mocpy.coverage_2d_from_ascii_str(value)
            return cls(cls.__create_key, index)
        elif format == "json":
            index = mocpy.coverage_2d_from_json_str(value)
            return cls(cls.__create_key, index)
        else:
            formats = ("ascii", "json")
            raise ValueError("format should be one of %s" % (str(formats)))
