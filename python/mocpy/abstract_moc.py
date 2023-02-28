
import os
from io import BytesIO
import numpy as np

from . import serializer

from . import mocpy

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


class AbstractMOC(serializer.IO):
    """Basic functions for manipulating MOCs."""

    _create_key = object()

    def __init__(self, create_key, create_sub_key, store_index):
        """Is an abstract coverage MOC.

        Args:
            create_key: Object ensuring __init__ is called by (sub)class-methods only
            create_sub_key: Object ensuring sub-classes __init__ is called from this (abstract) class or the sub-class
            store_index: index of the ST-MOC in the rust-side storage
        """
        assert (
            create_key == AbstractMOC._create_key
        ), "Abstract MOC instantiation is only allowed by sub-class methods"
        self.__create_sub_key = create_sub_key
        self._store_index = store_index

    def __repr__(self):
        return self.to_string(format="ascii", fold=80)

    def __del__(self):
        """Erase MOC."""
        mocpy.drop(self._store_index)

    def __eq__(self, other):
        """
        Test equality between this a MOC instance and ``another_moc`` of same type.

        Parameters
        ----------
        another_moc : `~mocpy.moc.MOC`, `~mocpy.tmoc.TMOC`, `~mocpy.stmoc.STMOC`
            The moc object to test the equality with

        Returns
        -------
        result : bool
            True if the self and ``another_moc`` are equal.
        """
        if not isinstance(other, AbstractMOC):
            raise TypeError(
                f"Cannot compare an AbstractMOC with a {type(other)}",
            )

        return mocpy.check_eq(self._store_index, other._store_index)

    def __copy__(self):
        mocpy.copy(self._store_index)
        return self.__class__(self.__create_sub_key, self._store_index)

    def __deepcopy__(self, memo):
        return self.__copy__()

    @staticmethod
    def store_index_dtype():
        usize_n_bits = mocpy.usize_n_bits()
        if usize_n_bits == 64:
            return np.uint64
        elif usize_n_bits == 32:
            return np.uint32
        else:
            raise ValueError("Unsupported store index usize type!")

    def __add__(self, moc):
        """
        Operator + definition.

        Computes the union of self with another MOC

        Parameters
        ----------
        moc : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            Another MOC to compute the union with.

        Returns
        -------
        result : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            The union of self and moc.
        """
        return self.union(moc)

    def __or__(self, moc):
        """
        Operator | definition.

        Computes the union of self with another MOC

        Parameters
        ----------
        moc : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            Another MOC to compute the union with.

        Returns
        -------
        result : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            The union of self and moc.
        """
        return self.union(moc)

    def __sub__(self, moc):
        """
        Operator - definition.

        Computes the difference of self with another MOC

        Parameters
        ----------
        moc : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            Another MOC to compute the difference with.

        Returns
        -------
        result : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            The difference of self with moc.
        """
        return self.difference(moc)

    def __and__(self, moc):
        """
        Operator & definition.

        Computes the intersection of self with another MOC

        Parameters
        ----------
        moc : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            Another MOC to compute the intersection with.

        Returns
        -------
        result : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            The intersection of self and moc.
        """
        return self.intersection(moc)

    def __invert__(self):
        """
        Unary operator ~ definition.

        Computes the complement of self

        Returns
        -------
        result : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            The complement MOC of self.
        """
        return self.complement()

    def empty(self):
        """
        Checks whether the MOC is empty or not.

        Returns
        -------
        result: bool
            True if the MOC instance is empty.
        """
        return mocpy.is_empty(self._store_index)

    @property
    def max_order(self):
        """Depth of the MOC instance."""
        raise NotImplementedError("Method max_order not implemented")

    @property
    def min_index(self):
        """Returns the smallest index (at the deepest absolute resolution) the MOC contains."""
        return mocpy.first_index(self._store_index)

    @property
    def max_index(self):
        """Returns the largest index (at the deepest absolute resolution) the MOC contains."""
        return mocpy.last_index(self._store_index)

    @property
    def uniq_gen(self):
        """
        Return a `np.array` of the generic uniq indices of the cell in the MOC.

        Warnings
        --------
        This is not defined in the MOC standard and is not HEALPix scpecific.

        Notes
        -----
        * It consists on the regular index with a sentinel bit placed at the immediate left
          of the index's MSB. At a given depth, the sentinel bit is always put o the same bit.
        * Because the uniq HEALPix encoding is not adapted for non-HEALPIx indices.
        """
        mocpy.to_uniq_gen(self._store_index)

    @property
    def uniq_zorder(self):
        """
        Return a `np.array` of the zorder uniq indices of the cell in the MOC.

        Warnings
        --------
        This is not defined in the MOC standard and is not HEALPix specific.

        Notes
        -----
        * It consists on a regular index shifted on the left so that indices at all level have the same MSB.
          Plus a sentinel bit placed at the immediate right of the LSB.
        * Because the uniq HEALPix encoding is not adapted for non-HEALPIx indices
          AND because the natural ordering of such indices follow the same order as the zorder indices
          (which is very useful for streaming processing, e.g. when dealing with multi-order maps)
        """
        mocpy.to_uniq_zorder(self._store_index)

    def flatten(self):
        """Return the list of indices of all cells in the MOC at the MOC depth."""
        return mocpy.flatten_to_moc_depth(self._store_index)

    def complement(self):
        """
        Returns the complement of the MOC instance.

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The resulting MOC.
        """
        index = mocpy.complement(self._store_index)

        return self.__class__(self.__create_sub_key, index)

    def intersection(self, another_moc, *args):
        """
        Intersection between the MOC instance and other MOCs.

        Parameters
        ----------
        another_moc : `~mocpy.moc.MOC`
            The MOC used for performing the intersection with self.
        args : `~mocpy.moc.MOC`
            Other additional MOCs to perform the intersection with.

        Returns
        -------
        result : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            The resulting MOC.
        """
        if args:
            store_indices = np.append(
                [self._store_index, another_moc._store_index],
                np.fromiter(
                    (arg._store_index for arg in args),
                    dtype=AbstractMOC.store_index_dtype(),
                ),
            )
            index = mocpy.multi_intersection(store_indices)
        else:
            index = mocpy.intersection(self._store_index, another_moc._store_index)

        return self.__class__(self.__create_sub_key, index)

    def union(self, another_moc, *args):
        """
        Union between the MOC instance and other MOCs.

        Parameters
        ----------
        another_moc : `~mocpy.moc.MOC`
            The MOC used for performing the union with self.
        args : `~mocpy.moc.MOC`
            Other additional MOCs to perform the union with.

        Returns
        -------
        result : `~mocpy.moc.MOC`/`~mocpy.tmoc.TimeMOC`
            The resulting MOC.
        """
        if args:
            store_indices = np.append(
                [self._store_index, another_moc._store_index],
                np.fromiter(
                    (arg._store_index for arg in args),
                    dtype=AbstractMOC.store_index_dtype(),
                ),
            )
            index = mocpy.multi_union(store_indices)
        else:
            index = mocpy.union(self._store_index, another_moc._store_index)

        return self.__class__(self.__create_sub_key, index)

    def symmetric_difference(self, another_moc, *args):
        """
        Symmetruic difference (XOR) between the MOC instance and other MOCs.

        Parameters
        ----------
        another_moc : `~mocpy.moc.MOC`
            The MOC used that will be substracted to self.
        args : `~mocpy.moc.MOC`
            Other additional MOCs to perform the difference with.

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The resulting MOC.
        """
        if args:
            store_indices = np.append(
                [self._store_index, another_moc._store_index],
                np.fromiter(
                    (arg._store_index for arg in args),
                    dtype=AbstractMOC.store_index_dtype(),
                ),
            )
            index = mocpy.multi_symmetric_difference(store_indices)
        else:
            index = mocpy.symmetric_difference(
                self._store_index, another_moc._store_index,
            )

        return self.__class__(self.__create_sub_key, index)

    def difference(self, another_moc, *args):
        """
        Difference between the MOC instance and other MOCs.

        Parameters
        ----------
        another_moc : `~mocpy.moc.MOC`
            The MOC used that will be substracted to self.
        args : `~mocpy.moc.MOC`
            Other additional MOCs to perform the difference with.

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The resulting MOC.
        """
        if args:
            store_indices = np.append(
                [another_moc._store_index],
                np.fromiter(
                    (arg._store_index for arg in args),
                    dtype=AbstractMOC.store_index_dtype(),
                ),
            )
            index_union = mocpy.multi_union(store_indices)
            index = mocpy.difference(self._store_index, index_union)
            mocpy.drop(index_union)
        else:
            index = mocpy.difference(self._store_index, another_moc._store_index)

        return self.__class__(self.__create_sub_key, index)

    def extended(self):
        """
        Returns the MOC extended by the external border made of cells at the MOC maximum depth.

        The only difference with respect to `add_neighbours` is that `extended` returns a new MOC
        instead of modifying the existing one.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The extended MOC
        """
        index = mocpy.extend(self._store_index)
        return self.__class__(self.__create_sub_key, index)

    def contracted(self):
        """
        Returns the MOC contracted by removing the internal border made of cells at the MOC maximum depth.

        The only difference with respect to `remove_neighbours` is that `contracted` returns a new MOC
        instead of modifying the existing one.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The extended MOC
        """
        index = mocpy.contract(self._store_index)
        return self.__class__(self.__create_sub_key, index)

    def add_neighbours(self):
        """
        Extends the MOC instance so that it includes the HEALPix cells touching its border.

        The depth of the HEALPix cells added at the border is equal to the maximum depth of the MOC instance.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            self extended by one degree of neighbours.
        """
        prev_store_index = self._store_index
        self._store_index = mocpy.extend(prev_store_index)
        # print(
        #    "\n@@ Manually change"
        #    + str(prev_store_index)
        #    + " into "
        #    + str(self._store_index)
        # )
        # print("\n@@ Manually drop" + str(prev_store_index))
        # mocpy.drop(prev_store_index)
        return self

    def remove_neighbours(self):
        """
        Removes from the MOC instance the HEALPix cells located at its border.

        The depth of the HEALPix cells removed is equal to the maximum depth of the MOC instance.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            self minus its HEALPix cells located at its border.
        """
        prev_store_index = self._store_index
        self._store_index = mocpy.contract(prev_store_index)
        # print(
        #    "\n@@ Manually change"
        #    + str(prev_store_index)
        #    + " into "
        #    + str(self._store_index)
        # )
        # print("\n@@ Manually drop" + str(prev_store_index))
        # mocpy.drop(prev_store_index)
        return self

    @classmethod
    def from_json(cls, json_moc):
        """
        Creates a MOC from a dictionary of HEALPix cell arrays indexed by their depth.

        Parameters
        ----------
        json_moc : dict(str : [int]
            A dictionary of HEALPix cell arrays indexed by their depth.

        Returns
        -------
        moc : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            the MOC.
        """
        import json

        return cls.from_string(
            json.dumps(json_moc, sort_keys=True, indent=2), format="json",
        )

    @classmethod
    def from_fits(cls, path_or_url):
        """
        Loads a MOC from a FITS file.

        The specified FITS file must store the MOC
        (i.e. the list of HEALPix cells it contains)
        in a binary HDU table.

        Parameters
        ----------
        path : str
            The path to the FITS file.

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The resulting MOC.

        """
        if isinstance(path_or_url, BytesIO):
            bytes = path_or_url.read()
            return cls._from_fits_raw_bytes(bytes)
        elif os.path.isfile(path_or_url):
            return cls.load(path_or_url, format="fits")
        else:
            import requests

            response = requests.get(path_or_url, headers={"User-Agent": "MOCPy"})
            if response:
                raw_bytes = BytesIO()
                raw_bytes.write(response.content)
                raw_bytes.seek(0)
                return cls.from_fits(raw_bytes)
            else:
                response.raise_for_status()

    @classmethod
    def from_str(cls, value):
        """
        Create a MOC from a str.

        This grammar is expressed is the `MOC IVOA <http://ivoa.net/documents/MOC/20190215/WD-MOC-1.1-20190215.pdf>`__
        specification at section 2.3.2.

        Parameters
        ----------
        value : str
            The MOC as a string following the grammar rules.

        Returns
        -------
        moc : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The resulting MOC

        Examples
        --------
        >>> from mocpy import MOC
        >>> moc = MOC.from_str("2/2-25 28 29 4/0 6/")
        """
        return cls.from_string(value, format="ascii")

    def degrade_to_order(self, new_order):
        """
        Degrades the MOC instance to a new, less precise, MOC.

        The maximum depth (i.e. the depth of the smallest cells that can be found in the MOC) of the
        degraded MOC is set to ``new_order``.

        Parameters
        ----------
        new_order : int
            Maximum depth of the output degraded MOC.

        Returns
        -------
        moc : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The degraded MOC.
        """
        raise NotImplementedError("Method degrade_to_order not implemented")

    def to_string(self, format="ascii", fold=0):
        """
        Writes the MOC into a string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        fold: int
            if >0, print ascii or json output with a maximum line width
        """
        if format == "ascii":
            if fold > 0:
                return mocpy.to_ascii_str_with_fold(self._store_index, fold)
            else:
                return mocpy.to_ascii_str(self._store_index)
        elif format == "json":
            if fold > 0:
                return mocpy.to_json_str_with_fold(self._store_index, fold)
            else:
                return mocpy.to_json_str(self._store_index)
        else:
            formats = ("ascii", "json")
            raise ValueError("format should be one of %s" % (str(formats)))

    def save(self, path, format="fits", overwrite=False, pre_v2=False, fold=0):
        """
        Writes the MOC to a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str or pathlib.Path
            The path to the file to save the MOC in.
        format : str, optional
            The format in which the MOC is saved.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        overwrite : bool, optional
            If the file already exists and you want to overwrite it, then set the  ``overwrite`` keyword.
            Default to False.
        fold: int
            if >0, print ascii or json output with a maximum line width
        """
        path = str(path)
        import os

        file_exists = os.path.isfile(path)

        if file_exists and not overwrite:
            raise OSError(
                "File {} already exists! Set ``overwrite`` to "
                "True if you want to replace it.".format(path),
            )

        if format == "fits":
            mocpy.to_fits_file(self._store_index, path, pre_v2)
        elif format == "ascii":
            if fold > 0:
                mocpy.to_ascii_file(self._store_index, path)
            else:
                mocpy.to_ascii_file_with_fold(self._store_index, path, fold)
        elif format == "json":
            if fold > 0:
                mocpy.to_json_file(self._store_index, path)
            else:
                mocpy.to_json_file_with_fold(self._store_index, path, fold)
        else:
            formats = ("fits", "ascii", "json")
            raise ValueError("format should be one of %s" % (str(formats)))
