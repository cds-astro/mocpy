# -*- coding: utf-8 -*
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from astropy.io import fits
from astropy.table import Table

from .interval_set import IntervalSet
from . import utils, serializer

from . import mocpy

__author__ = "Thomas Boch, Matthieu Baumann, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, baumannmatthieu0@gmail.com, francois-xavier.pineau@astro.unistra.fr"


class AbstractMOC(serializer.IO):
    """
    Basic functions for manipulating MOCs.
    """
    
    
    __LARK_PARSER_STR = None

    def __init__(self, interval_set=None):
        interval = IntervalSet() if interval_set is None else interval_set
        self._interval_set = interval
        self._fits_column_name = 'UNIQ'

    def __repr__(self):
        return self._interval_set.__repr__()

    def __eq__(self, another_moc):
        """
        Test equality between thr MOC instance and ``another_moc``

        Parameters
        ----------
        another_moc : `~mocpy.moc.MOC`
            The moc object to test the equality with

        Returns
        -------
        result : bool
            True if the interval sets of self and ``another_moc`` are equal (the interval sets are checked
            for consistency before comparing them).
        """
        if not isinstance(another_moc, AbstractMOC):
            raise TypeError('Cannot compare an AbstractMOC with a {0}'.format(type(another_moc)))

        return self._interval_set == another_moc._interval_set

    def empty(self):
        """
        Checks whether the MOC is empty.

        A MOC is empty when its list of HEALPix cell ranges is empty.

        Returns
        -------
        result: bool
            True if the MOC instance is empty.
        """
        return self._interval_set.empty()

    @property
    def max_order(self):
        """
        Depth of the smallest cells found in the MOC instance.
        """
        #depth = mocpy.coverage_depth(self._interval_set._intervals)
        #depth = np.uint8(depth)
        #return depth
        raise NotImplementedError("Method max_order not implemented")


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
        intervals = self._interval_set.intersection(another_moc._interval_set)

        for moc in args:
            intervals = intervals.intersection(moc._interval_set)

        return self.__class__(intervals)

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
        intervals = self._interval_set.union(another_moc._interval_set)

        for moc in args:
            intervals = intervals.union(moc._interval_set)

        return self.__class__(intervals)

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
        intervals = self._interval_set.difference(another_moc._interval_set)

        for moc in args:
            intervals = intervals.difference(moc._interval_set)

        return self.__class__(intervals)

    def complement(self):
        """
        Returns the complement of the MOC instance.

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The resulting MOC.
        """
        # intervals = self._interval_set.complement()
        # return self.__class__(intervals)
        raise NotImplementedError("Method complement not implemented")

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
        intervals = mocpy.coverage_from_json(json_moc)
        return cls(IntervalSet(intervals, make_consistent=False))

    @classmethod
    def from_fits(cls, filename):
        """
        Loads a MOC from a FITS file.

        The specified FITS file must store the MOC
        (i.e. the list of HEALPix cells it contains)
        in a binary HDU table.

        Parameters
        ----------
        filename : str
            The path to the FITS file.

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            The resulting MOC.
            
        Warning
        -------
        Not compatible with the last version of the MOC 2.0 standard.
            
        """
        import warnings
        warnings.warn('This method is deprecated. Use MOC.load(path, "fits") instead!', DeprecationWarning)
        table = Table.read(filename)
        first_column_index = table.colnames[0]
        intervals = mocpy.to_nested(table[first_column_index].astype(np.uint64))
        return cls(IntervalSet(intervals, make_consistent=False))

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
        import warnings
        warnings.warn('This method is deprecated. Use MOC.load(path, "ascii") instead!', DeprecationWarning)
        
        
        # Import lark parser when from_str is called
        # at least one time
        from lark import Lark, Transformer
        class ParsingException(Exception):
            pass

        class TreeToJson(Transformer):
            def value(self, items):
                res = {}
                for item in items:
                    if item is not None: # Do not take into account the "sep" branches
                        res.update(item)
                return res

            def sep(self, items):
                pass

            def depthpix(self, items):
                depth = str(items[0])
                pixs_l = items[1:][0]
                return {depth: pixs_l}

            def uniq_pix(self, pix):
                if pix:
                    return [int(pix[0])]

            def range_pix(self, range_pix):
                lower_bound = int(range_pix[0])
                upper_bound = int(range_pix[1])
                return np.arange(start=lower_bound, stop=upper_bound + 1, dtype=int).tolist()

            def pixs(self, items):
                ipixs = []
                for pix in items:
                    if pix is not None: # Do not take into account the "sep" branches
                        ipixs.extend(pix)
                return ipixs

        # Initialize the parser when from_str is called
        # for the first time
        if AbstractMOC.__LARK_PARSER_STR is None:
            AbstractMOC.__LARK_PARSER_STR = Lark(r"""
                value: depthpix (sep+ depthpix)*
                depthpix : INT "/" sep* pixs
                pixs : pix (sep+ pix)*
                pix : INT? -> uniq_pix
                    | (INT "-" INT) -> range_pix
                sep : " " | "\n" | "\r"

                %import common.INT
                """, start='value')

        try:
            tree = AbstractMOC.__LARK_PARSER_STR.parse(value)
        except Exception as err:
            raise ParsingException("Could not parse {0}. \n Check the grammar \
                section 2.3.2 of http://ivoa.net/documents/MOC/20190215/WD-MOC-1.1-20190215.pdf \
                to see the correct syntax for writing a MOC from a str".format(value))

        moc_json = TreeToJson().transform(tree)
        return cls.from_json(moc_json)

    def _uniq_format(self):
        return self._interval_set.uniq

    @staticmethod
    def _to_json(uniq):
        """
        Serializes a MOC to the JSON format.

        Parameters
        ----------
        uniq : `~numpy.ndarray`
            The array of HEALPix cells representing the MOC to serialize.

        Returns
        -------
        result_json : {str : [int]}
            A dictionary of HEALPix cell lists indexed by their depth.
        """
        result_json = {}

        depth, ipix = utils.uniq2orderipix(uniq)
        min_depth = np.min(depth[0])
        max_depth = np.max(depth[-1])

        for d in range(min_depth, max_depth + 1):
            pix_index = np.where(depth == d)[0]
            if pix_index.size:
                # there are pixels belonging to the current order
                ipix_depth = ipix[pix_index]
                result_json[str(d)] = ipix_depth.tolist()

        return result_json
        #return mocpy.coverage_to_json(intervals)

    @staticmethod
    def _to_str(uniq):
        """
        Serializes a MOC to the STRING format.

        HEALPix cells are separated by a comma. The HEALPix cell at order 0 and number 10 is encoded
        by the string: "0/10", the first digit representing the depth and the second the HEALPix cell number
        for this depth. HEALPix cells next to each other within a specific depth can be expressed as a range and 
        therefore written like that: "12/10-150". This encodes the list of HEALPix cells from 10 to 150 at the
        depth 12.

        Parameters
        ----------
        uniq : `~numpy.ndarray`
            The array of HEALPix cells representing the MOC to serialize.

        Returns
        -------
        result : str
            The serialized MOC.
        """
        def write_cells(serial, a, b, sep=''):
            if a == b:
                serial += '{0}{1}'.format(a, sep)
            else:
                serial += '{0}-{1}{2}'.format(a, b, sep)
            return serial

        res = ''

        if uniq.size == 0: 
            return res

        depth, ipixels = utils.uniq2orderipix(uniq)
        min_depth = np.min(depth[0])
        max_depth = np.max(depth[-1])

        for d in range(min_depth, max_depth + 1):
            pix_index = np.where(depth == d)[0]

            if pix_index.size > 0:
                # Serialize the depth followed by a slash
                res += '{0}/'.format(d)

                # Retrieve the pixel(s) for this depth
                ipix_depth = ipixels[pix_index]
                if ipix_depth.size == 1:
                    # If there is only one pixel we serialize it and
                    # go to the next depth
                    res = write_cells(res, ipix_depth[0], ipix_depth[0])
                else:
                    # Sort them in case there are several
                    ipix_depth = np.sort(ipix_depth)

                    beg_range = ipix_depth[0]
                    last_range = beg_range

                    # Loop over the sorted pixels by tracking the lower bound of
                    # the current range and the last pixel.
                    for ipix in ipix_depth[1:]:
                        # If the current pixel does not follow the previous one
                        # then we can end a range and serializes it
                        if ipix > last_range + np.uint64(1):
                            res = write_cells(res, beg_range, last_range, sep=' ')
                            # The current pixel is the beginning of a new range
                            beg_range = ipix

                        last_range = ipix

                    # Write the last range
                    res = write_cells(res, beg_range, last_range)

                # Add a ' ' separator before writing serializing the pixels of the next depth
                res += ' '

        # Remove the last ' ' character
        res = res[:-1]

        return res

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
        #intervals = mocpy.coverage_degrade(self._interval_set._intervals, new_order)
        #return self.__class__(IntervalSet(intervals, make_consistent=False))
        raise NotImplementedError("Method degrade_to_order not implemented")

    def refine_to_order(self, min_depth):
        """
        Internal method
        """
        #intervals = mocpy.coverage_merge_nested_intervals(self._interval_set._intervals, min_depth)
        #interval_set = IntervalSet(intervals, make_consistent=False)
        #return self.__class__(interval_set)
        raise NotImplementedError("Method refine_to_order not implemented")

