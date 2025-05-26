import functools

import numpy as np
from astropy import units as u

from .. import MOC, FrequencyMOC, mocpy
from ..abstract_moc import AbstractMOC

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = (
    "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, "
    "manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"
)


def validate_frequencies_and_positions(function):
    """Validate that the entries are either scalar of 1D array of values.

    Parameters
    ----------
    function : <class 'function'>
        must have the signature function(self, my_thing_to_test, **kwargs)
    """

    @functools.wraps(function)
    def _validate_frequencies_and_positions_wrap(self, *args, **kwargs):
        # convert to at least 1D
        args = list(map(np.atleast_1d, args))
        kwargs = {
            key: (np.atleast_1d(value) if "order" not in key else value)
            for key, value in kwargs.items()
        }
        from_keywords_arguments = [
            value for key, value in kwargs.items() if "order" not in key
        ]
        # test ndim
        if any(arg.ndim != 1 for arg in args + from_keywords_arguments):
            raise ValueError(
                "Frequencies and positions must be scalar quantities or 1D arrays."
            )
        # test that they all have the same length
        if len({len(arg) for arg in args + from_keywords_arguments}) != 1:
            raise ValueError("Frequencies and positions must have the same lengths.")
        return function(self, *args, **kwargs)

    return _validate_frequencies_and_positions_wrap


class SFMOC(AbstractMOC):
    """Spatial-Frequency Coverage class."""

    def __init__(self, store_index):
        """Is a Spatial-Frequency Coverage (SF-MOC).

        Args:
            store_index: index of the SF-MOC in the rust-side storage
        """
        self.store_index = store_index

    def __del__(self):
        """Erase SFMOC."""
        super().__del__()

    def __eq__(self, other):
        """Assert equality between SF-MOCs."""
        return super().__eq__(other)

    @property
    def max_order(self):
        """Maximum order of the SF-MOC.

        Returns
        -------
        (int, int)
            (max_order_frequency, max_order_space)
        """
        return mocpy.coverage_sf_depth(self.store_index)

    @property
    def min_frequency(self):
        """Return SFMOC min frequency."""
        return mocpy.coverage_sf_min_freq(self.store_index) * u.Hz

    @property
    def max_frequency(self):
        """Return SFMOC max frequency."""
        return mocpy.coverage_sf_max_freq(self.store_index) * u.Hz

    @classmethod
    def n_cells(cls, order, *, dimension):
        """Get the number of cells for a given order.

        Parameters
        ----------
        order : int
            The order. It is comprised between 0 and `~mocpy.moc.MOC.MAX_ORDER` if
            dimension='space' and between 0 and `~mocpy.tmoc.FrequencyMOC.MAX_ORDER` if
            dimension='frequency'.

        dimension : str
            Can be either 'frequency' or 'space'.

        Returns
        -------
        int
            The number of cells at the given order

        Examples
        --------
        >>> from mocpy import SFMOC
        >>> SFMOC.n_cells(0, dimension='space')
        12
        """
        if dimension == "space":
            return MOC.n_cells(order)
        if dimension == "frequency":
            return FrequencyMOC.n_cells(order)
        raise ValueError(
            f"Dimension should be either 'time' of 'frequency' but '{dimension}' was provided.",
        )

    @classmethod
    def new_empty(cls, max_order_frequency, max_order_space):
        """Create a new empty SFMOC.

        Parameters
        ----------
        max_order_frequency : int
            The frequency resolution of the SFMOC. Should be comprised between 0 and 59.
        max_order_space : int
            The space resolution of the SFMOC. Should be comprised between 0 and 29.

        Returns
        -------
        `~mocpy.SFMOC`

        Examples
        --------
        >>> from mocpy import SFMOC
        >>> SFMOC.new_empty(20, 12)
        f20/ s12/

        """
        index = mocpy.new_empty_sfmoc(
            np.uint8(max_order_frequency),
            np.uint8(max_order_space),
        )
        return cls(index)

    def is_empty(self):
        """Check whether the Space-Frequency coverage is empty.

        Returns
        -------
        bool

        Examples
        --------
        >>> from mocpy import SFMOC
        >>> sfmoc = SFMOC.new_empty(20, 12)
        >>> sfmoc.is_empty()
        True
        """
        return mocpy.is_empty(self.store_index)

    @classmethod
    @validate_frequencies_and_positions
    def from_frequencies_and_positions(
        cls, frequencies, lon, lat, *, max_order_frequency, max_order_space
    ):
        """
        Create a Space-Frequency Coverage from a set of frequencies and positions.

        Parameters
        ----------
        frequencies : `astropy.units.Quantity`
            An astropy Quantity of physical type 'frequency'
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates corresponding to the frequencies.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates corresponding to the frequencies.
        max_order_frequency : int
            Frequency order. Should be comprised between 0 and 59.
        max_order_space : int
            Spatial order.

        Returns
        -------
        `~mocpy.SFMOC`

        Examples
        --------
        >>> import astropy.units as u
        >>> from mocpy import SFMOC
        >>> frequencies = [1, 2, 3] * u.Hz
        >>> lon = [0, 1, 2] * u.deg
        >>> lat = [0, 1, 2] * u.deg
        >>> sfmoc = SFMOC.from_frequencies_and_positions(frequencies, lon, lat,
        ...                                              max_order_frequency=20,
        ...                                              max_order_space=12)
        >>> sfmoc
        f20/770048
        s12/79691776
        f20/778240
        s12/79697029
        f20/782336
        s12/79712788
        f20/ s12/
        """
        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)
        frequencies = frequencies.to_value("Hz").astype(np.float64)

        return cls(
            mocpy.from_freq_lonlat(
                frequencies, max_order_frequency, lon, lat, max_order_space
            )
        )

    @classmethod
    @validate_frequencies_and_positions
    def from_frequency_ranges_and_positions(
        cls,
        frequencies_min,
        frequencies_max,
        lon,
        lat,
        *,
        max_order_frequency,
        max_order_space,
    ):
        """Create a SF coverage from a range of frequencies for each position.

        Parameters
        ----------
        frequencies_min : `astropy.units.Quantity`
            An astropy Quantity of physical type 'frequency'
        frequencies_max : `astropy.units.Quantity`
            An astropy Quantity of physical type 'frequency'
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed at a specific time.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed at a specific time.
        max_order_frequency : int
            Frequency order.
        max_order_space : int
            Spatial order.

        Returns
        -------
        `~mocpy.SFMOC`
            The resulting Space-Frequency Coverage map.

        Examples
        --------
        >>> import astropy.units as u
        >>> from mocpy import SFMOC
        >>> frequencies_min = [0.01, 0.02, 0.03] * u.Hz
        >>> frequencies_max = [0.1, 0.2, 0.3] * u.Hz
        >>> lon = [0, 1, 2] * u.deg
        >>> lat = [0, 1, 2] * u.deg
        >>> sfmoc = SFMOC.from_frequency_ranges_and_positions(frequencies_min,
        ...                                                   frequencies_max, lon, lat,
        ...                                                   max_order_frequency=10,
        ...                                                   max_order_space=12)
        >>> sfmoc
        f8/175
        9/349 352
        10/
        s12/79691776
        f9/353-354
        10/710
        s12/79691776 79697029
        f7/89
        8/180
        10/711 724
        s12/79691776 79697029 79712788
        f8/182
        9/363
        10/725 732
        s12/79697029 79712788
        f9/367-368
        10/733
        s12/79712788
        f10/ s12/
        """
        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)
        frequencies_min = frequencies_min.to_value("Hz").astype(np.float64)
        frequencies_max = frequencies_max.to_value("Hz").astype(np.float64)

        index = mocpy.from_freq_ranges_lonlat(
            frequencies_min,
            frequencies_max,
            max_order_frequency,
            lon,
            lat,
            max_order_space,
        )

        return cls(index)

    @classmethod
    @validate_frequencies_and_positions
    def from_spatial_coverages(
        cls,
        frequencies_min,
        frequencies_max,
        spatial_coverages,
        *,
        max_order_frequency,
    ):
        """Create a ST coverage from frequency ranges associated to spatial coverages.

        Parameters
        ----------
        frequencies_min : `~astropy.units.Quantity`
            An astropy Quantity of physical type 'frequency'
        frequencies_max : `~astropy.units.Quantity`
            An astropy Quantity of physical type 'frequency'
        spatial_coverages : list[`~mocpy.MOC`]
            List of spatial coverages.
        max_order_frequency : int
            Frequency order.

        Returns
        -------
        `~mocpy.SFMOC`

        Examples
        --------
        >>> import astropy.units as u
        >>> from mocpy import MOC, SFMOC
        >>> sfmoc = SFMOC.from_spatial_coverages(
        ...             frequencies_min=[0.1]*u.Hz, frequencies_max=[10]*u.Hz,
        ...             spatial_coverages=[MOC.from_string("5/14-21")],
        ...             max_order_frequency=10)
        >>> sfmoc
        f5/23
        7/91 96
        8/181
        9/388
        10/
        s4/4
        5/14-15 20-21
        f10/ s5/
        """
        spatial_coverages_indices = np.fromiter(
            (arg.store_index for arg in spatial_coverages),
            dtype=AbstractMOC._store_index_dtype(),
        )
        index = mocpy.from_frequency_ranges_spatial_coverages(
            frequencies_min,
            frequencies_max,
            max_order_frequency,
            spatial_coverages_indices,
        )

        return cls(index)

    def query_by_frequency(self, fmoc):
        """Query the SF-MOC by frequency F-MOC.

        This will perform the union of all the spatial coverages lying in a set of time ranges.

        Parameters
        ----------
        fmoc : ~mocpy.FrequencyMOC``
            Frequency MOC.

        Returns
        -------
        `~mocpy.MOC`
            The spatial coverage being observed within the input frequency ranges

        Examples
        --------
        >>> from mocpy import MOC, FrequencyMOC, SFMOC
        >>> sfmoc = SFMOC.from_string('''
        ... f15/0-10
        ... s12/0-100
        ... f15/11-20
        ... s12/101-200
        ... ''')
        >>> fmoc = FrequencyMOC.from_string("15/0-2")
        >>> MOC.from_sfmoc_frequency_fold(fmoc, sfmoc)
        9/0
        10/4-5
        11/24
        12/100
        """
        return MOC.from_sfmoc_frequency_fold(fmoc, self)

    def query_by_space(self, smoc):
        """Query the SF-MOC by space coverage.

        This will perform the union of all the frequency ranges which associated
        spatial coverages fall within the given spatial MOC.

        Parameters
        ----------
        smoc : `~mocpy.MOC`
            The spatial coverage.

        Returns
        -------
        `~mocpy.FrequencyMOC`
            The Frequency coverage corresponding to the Spatial MOC

        Examples
        --------
        >>> from mocpy import MOC, SFMOC, FrequencyMOC as FMOC
        >>> sfmoc = SFMOC.from_string('''f10/0-20
        ... s12/0-100
        ... f10/21-40
        ... s12/101-200
        ... ''')
        >>> moc = MOC.from_string("12/0-100")
        >>> fmoc = FMOC.from_sfmoc_space_fold(sfmoc, moc)
        >>> fmoc
        6/0
        8/4
        10/20
        """
        return FrequencyMOC.from_sfmoc_space_fold(self, smoc)

    @validate_frequencies_and_positions
    def contains(self, frequencies, lon, lat):
        """Test if the Frequency-Space combinations fall within the SFMOC.

        Parameters
        ----------
        frequencies : `~astropy.units.Quantity`
            Astropy quantities of physical type ``frequency``.
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed with a specific frequency.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed with a specific frequency.

        Returns
        -------
        `~np.array`
            A boolean array with the same length than the given Frequency-Position
            couples. True if the Frequency-Position is within the SFMOC, False
            otherwise.

        Examples
        --------
        >>> from mocpy import SFMOC, MOC
        >>> import astropy.units as u
        >>> moc = MOC.from_cone(0*u.deg, 0*u.deg, radius=10*u.deg, max_depth=10)
        >>> sfmoc = SFMOC.from_spatial_coverages(0.01*u.Hz, 100*u.Hz,
        ...                                      moc, max_order_frequency=40)
        >>> # one inside, one outside
        >>> sfmoc.contains([10, 10000]*u.Hz, [0.1, 20]*u.deg, [0.1, 20]*u.deg)
        array([ True, False])
        """
        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)
        frequencies = frequencies.to_value("Hz").astype(np.float64)
        return mocpy.sfmoc_contains(self.store_index, frequencies, lon, lat)

    @classmethod
    # A002: Argument `format` is shadowing a python function
    def load(cls, path, format="fits"):  # noqa: A002
        """Load the Space-Frequency MOC from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not
        officially supported by the IVOA.

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
        `~mocpy.SFMOC`
        """
        path = str(path)
        if format == "fits":
            index = mocpy.coverage_sf_from_fits_file(path)
            return cls(index)
        if format == "ascii":
            index = mocpy.coverage_sf_from_ascii_file(path)
            return cls(index)
        if format == "json":
            index = mocpy.coverage_sf_from_json_file(path)
            return cls(index)
        formats = ("fits", "ascii", "json")
        raise ValueError(f"format should be one of {formats}")

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load MOC from raw bytes of a FITS file."""
        index = mocpy.sfmoc_from_fits_raw_bytes(raw_bytes)
        return cls(index)

    @classmethod
    # A002: Argument `format` is shadowing a python function
    def from_string(cls, value, format="ascii"):  # noqa: A002
        """
        Deserialize the Spatial MOC from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially
        supported by the IVOA.

        Parameters
        ----------
        format : str, optional
            The format in which the MOC is serialized.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".

        Returns
        -------
        `~mocpy.SFMOC`
            The ST-MOC build from the given string.
        """
        if format == "ascii":
            index = mocpy.coverage_sf_from_ascii_str(value)
            return cls(index)
        if format == "json":
            index = mocpy.coverage_sf_from_json_str(value)
            return cls(index)
        formats = ("ascii", "json")
        raise ValueError(f"format should be one of {formats}")
