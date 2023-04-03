import numpy as np
from astropy import units as u

from ..abstract_moc import AbstractMOC

from .. import mocpy

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


class FrequencyMOC(AbstractMOC):
    """Multi-order frequency coverage class. Experimental."""

    # Maximum order of F-MOCs
    MAX_ORDER = np.uint8(59)
    # Upper limit (exclusive) on a F-MOC index = 2^60
    MAX_INDEX_EXCLUSIVE = 1152921504606846976
    # Smallest value, in Hz, a F-MOC can encode
    FREQ_MIN_HZ = 5.048_709_793_414_476e-29
    # Largest value, in Hz, a F-MOC can encode
    FREQ_MAX_HZ = 5.846_006_549_323_611e48

    __create_key = object()

    def __init__(self, create_key, store_index):
        """Is a Frequency Coverage (S-MOC).

        Args:
            create_key: Object ensure __init__ is called by super-class/class-methods only
            store_index: index of the S-MOC in the rust-side storage
        """
        super().__init__(
            AbstractMOC._create_key,
            FrequencyMOC.__create_key,
            store_index,
        )
        if create_key != FrequencyMOC.__create_key:
            raise PermissionError(
                "F-MOC instantiation is only allowed by class or super-class methods",
            )

    @property
    def max_order(self):
        """Depth/order of the F-MOC."""
        depth = mocpy.get_fmoc_depth(self._store_index)
        return np.uint8(depth)

    def to_hz_ranges(self):
        """Return the Hertz ranges this FrequencyMOC contains, in Hertz."""
        return np.asarray(
            mocpy.to_freq_ranges(self._store_index) * u.Hz,
            dtype=np.float64,
        )

    @property
    def to_depth59_ranges(self):
        """Return the list of ranges this FrequencyMOC contains."""
        return mocpy.to_ranges(self._store_index)

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
        moc : `~mocpy.fmoc.FrequencyMOC`
            The degraded MOC.
        """
        index = mocpy.degrade(self._store_index, new_order)
        return FrequencyMOC(FrequencyMOC.__create_key, index)

    @classmethod
    def new_empty(cls, max_depth):
        """
        Create a new empty FrequencyMOC of given depth.

        Parameters
        ----------
        max_depth : int, The resolution of the FrequencyMOC


        Returns
        -------
        moc : `~mocpy.fmoc.FrequencyMOC`
            The MOC
        """
        index = mocpy.new_empty_fmoc(np.uint8(max_depth))
        return cls(cls.__create_key, index)

    @classmethod
    def from_depth59_ranges(cls, order, ranges):
        """
        Create a FrequencyMOC from a set of FrequencyMOC ranges at order 59.

        Parameters
        ----------
        order : int, The resolution of the FrequencyMOC
        ranges: `~numpy.ndarray`
                 a N x 2 numpy array representing the set of depth 61 ranges.

        Returns
        -------
        moc : `~mocpy.tmoc.TimeMOC`
            The MOC
        """
        ranges = np.zeros((0, 2), dtype=np.uint64) if ranges is None else ranges

        if ranges.shape[1] != 2:
            raise ValueError(
                f"Expected a N x 2 numpy ndarray but second dimension is {ranges.shape[1]}",
            )

        if ranges.dtype is not np.uint64:
            ranges = ranges.astype(np.uint64)

        index = mocpy.from_fmoc_ranges_array2(np.uint8(order), ranges)
        return cls(cls.__create_key, index)

    @classmethod
    def from_frequencies(cls, order, frequencies):
        """
        Create a FrequencyMOC from a `astropy.units.Quantity` that are internally converted in `Hz`.

        Parameters
        ----------
        order : int, The resolution of the FrequencyMOC: see `relative_precision_to_order`
        frequencies : `astropy.units.Quantity`
            Quantity converted internally in `Hz`

        Returns
        -------
        frequency_moc : `~mocpy.fmoc.FrequencyMOC`
        """
        frequencies = frequencies.to(u.Hz)
        frequencies = np.atleast_1d(frequencies)

        store_index = mocpy.from_freq_values(order, frequencies)
        return cls(cls.__create_key, store_index)

    @classmethod
    def from_frequency_ranges(cls, order, min_freq, max_freq):
        """
        Create a FrequencyMOC from a range defined by two `astropy.units.Quantity` converted in `Hz`.

        Parameters
        ----------
        order : int, The resolution of the FrequencyMOC: see `relative_precision_to_order`
        min_freq : `astropy.units.Quantity`
            astropy quantity converted in Hz and defining the left part of the intervals
        max_freq : `astropy.units.Quantity`
            astropy quantity converted in Hz and defining the right part of the intervals

        Returns
        -------
        frequency_moc : `~mocpy.fmoc.FrequencyMOC`
        """
        min_freq = min_freq.to(u.Hz)
        min_freq = np.atleast_1d(min_freq)

        max_freq = max_freq.to(u.Hz)
        max_freq = np.atleast_1d(max_freq)

        if min_freq.shape != max_freq.shape:
            raise ValueError(
                f"Mismatch between min_freq and max_freq of shapes {min_freq.shape} and {max_freq.shape}",
            )

        store_index = mocpy.from_freq_ranges(
            order,
            min_freq,
            max_freq,
        )
        return cls(cls.__create_key, store_index)

    @property
    def min_freq(self):
        """
        Get the `~astropy.units.Quantity` frequency of the fmoc first observation.

        Returns
        -------
        min_freq : `astropy.units.Quantity`
            frequency of the first observation

        """
        return np.atleast_1d(mocpy.first_fmoc_hz(self._store_index) * u.Hz)

    @property
    def max_freq(self):
        """
        Get the `~astropy.units.Quantity` frequency of the fmoc last observation.

        Returns
        -------
        max_freq : `~astropy.units.Quantity`
            frequency of the last observation

        """
        return np.atleast_1d(mocpy.last_fmoc_hz(self._store_index) * u.Hz)

    def contains(self, frequencies, keep_inside=True):
        """
        Get a mask array (e.g. a numpy boolean array) of times being inside (or outside) the FMOC instance.

        Parameters
        ----------
        frequencies : `astropy.units.Quantity`
            astropy quantity (converted into Hz) to check whether they are contained in the FMOC or not.
        keep_inside : bool, optional
            True by default. If so the filtered table contains only observations that are located the MOC.
            If ``keep_inside`` is False, the filtered table contains all observations lying outside the MOC.

        Returns
        -------
        array : `~numpy.darray`
            A mask boolean array
        """
        freq = frequencies.to(u.Hz)
        freq = np.atleast_1d(freq)

        mask = mocpy.filter_freq(self._store_index, freq)

        if keep_inside:
            return mask
        return ~mask

    @staticmethod
    def order_to_relative_precision(order):
        """
        Convert a Frequency Moc order to a relative precision range.

        Parameters
        ----------
        order : int
            order to convert

        Returns
        -------
        rel_prec : relative precision
            time equivalent to ``order``

        """
        return np.power([2.0, 2.0], [6 - order, 7 - order])

    @staticmethod
    def relative_precision_to_order(frequency_precision):
        """
        Convert a relative frequency precision into a FrequencyMOC order.

        Parameters
        ----------
        frequency_precision : `float`
            precision to be converted in Frequency MOC order

        Returns
        -------
        order : int
            The less precise order fulfilling the given ``frequency_precision``.

        """
        order = int(7 - np.log(frequency_precision) / np.log(2))
        if order > FrequencyMOC.MAX_ORDER:
            order = FrequencyMOC.MAX_ORDER
        elif order < 0:
            order = 0
        return np.uint8(order)

    @classmethod
    def load(cls, path, format="fits"):  # noqa: A002
        """
        Load the FrequencyMOC MOC from a file.

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
            index = mocpy.frequency_moc_from_fits_file(path)
            return cls(cls.__create_key, index)
        if format == "ascii":
            index = mocpy.frequency_moc_from_ascii_file(path)
            return cls(cls.__create_key, index)
        if format == "json":
            index = mocpy.frequency_moc_from_json_file(path)
            return cls(cls.__create_key, index)
        formats = ("fits", "ascii", "json")
        raise ValueError("format should be one of %s" % (str(formats)))

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load Frequency MOC from raw bytes of a FITS file."""
        index = mocpy.frequency_moc_from_fits_raw_bytes(raw_bytes)
        return cls(cls.__create_key, index)

    @classmethod
    def from_string(cls, value, format="ascii"):  # noqa: A002
        """
        Deserialize the Frequency MOC from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        WARNING: the serialization must be strict, i.e. **must not** contain overlapping elements

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == "ascii":
            index = mocpy.frequency_moc_from_ascii_str(value)
            return cls(cls.__create_key, index)
        if format == "json":
            index = mocpy.frequency_moc_from_json_str(value)
            return cls(cls.__create_key, index)
        formats = ("ascii", "json")
        raise ValueError("format should be one of %s" % (str(formats)))
