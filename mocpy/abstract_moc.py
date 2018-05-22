#!/usr/bin/env python
# -*- coding: utf-8 -*

"""
abstract_moc.py

basic functions for manipulating mocs

"""

from __future__ import absolute_import, division, print_function, unicode_literals
from . import py23_compat

import numpy as np
from astropy.io import fits
from astropy_healpix.healpy import ring2nest

from .interval_set import IntervalSet
from . import utils

__author__ = "Thomas Boch, Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, matthieu.baumann@astro.unistra.fr"


class AbstractMoc:
    HPY_MAX_NORDER = 29

    def __init__(self):
        self._interval_set = IntervalSet()

    def __repr__(self):
        return self._interval_set.__repr__()

    def __eq__(self, another_moc):
        """
        Test equality between self and ``another_moc``

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMoc`
            the moc object to test the equality with

        Returns
        -------
        result : bool
            True if the interval sets of self and ``another_moc`` are equal (the interval sets are checked
            for consistency before comparing them).

        """

        if not isinstance(another_moc, AbstractMoc):
            raise TypeError('Cannot compare an AbstractMoc with a {0}'.format(type(another_moc)))

        result = (self._interval_set.intervals == another_moc._interval_set.intervals)
        return result

    def add_pix_list(self, order, i_pix_l, nest=True):
        """
        Add a list of HEALPix pixels at a given order to the
        current mocpy object (This method uses numpy vectorization)

        Parameters
        ----------
        order : int
            The order of the pixel list to add
        i_pix_l : [int]
            The pixel list at ``order``
        nest : bool, optional
            True by default

        """

        add_pix_vect = np.vectorize(self.add_pix)
        add_pix_vect(order, i_pix_l, nest=nest)

    def add_pix(self, order, i_pix, nest=True):
        """
        Add a HEALPix pixel at a given order
        to the current mocpy object (This method uses numpy vectorization)

        Parameters
        ----------
        order : int
            The order of the pixel to add
        i_pix : int
            The pixel at ``order``
        nest : bool, optional
            True by default

        """

        if order > AbstractMoc.HPY_MAX_NORDER:
            raise ValueError('order can not be greater than AbstractMoc.HPY_MAX_NORDER')

        if not nest:
            i_pix = ring2nest(1 << order, i_pix)

        p1 = i_pix
        p2 = i_pix + 1
        shift = 2 * (AbstractMoc.HPY_MAX_NORDER - order)

        self._interval_set.add((p1 << shift, p2 << shift))

    @property
    def max_order(self):
        """
        This returns the deepest order needed to describe the current _interval_set

        """

        #  TODO: cache value
        combo = int(0)
        for iv in self._interval_set.intervals:
            combo |= iv[0] | iv[1]

        ret = AbstractMoc.HPY_MAX_NORDER - (utils.number_trailing_zeros(combo) // 2)
        if ret < 0:
            ret = 0

        return ret

    def _ensure_consistency(self):
        """
        Internal operation on the interval set of the mocpy object ensuring its consistency

        This private method is called when the user asks for :
        * plotting the moc
        * performing simple operations on it : intersection, union, difference of mocs
        * serializing the moc into json/fits format

        """
        self._interval_set.intervals

    def intersection(self, another_moc):
        """
        Intersection between self and moc

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMoc`
            the MOC/TimeMoc used for performing the intersection with self

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self & ``moc``

        """

        iv_set_intersection = self._interval_set.intersection(another_moc._interval_set)
        result = self.__class__.from_interval_set(iv_set_intersection)
        return result

    def union(self, another_moc):
        """
        Union between self and moc

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMoc`
            the MOC/TimeMoc to bind to self

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self | ``moc``

        """

        iv_set_union = self._interval_set.union(another_moc._interval_set)
        result = self.__class__.from_interval_set(iv_set_union)
        return result

    def difference(self, moc):
        """
        Difference between self and moc

        Parameters
        ----------
        moc : `~mocpy.abstract_moc.AbstractMoc`
            the MOC/TimeMoc to substract from self

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self - ``moc``

        """

        iv_set_difference = self._interval_set.difference(moc._interval_set)
        result = self.__class__.from_interval_set(iv_set_difference)
        return result

    def complement(self):
        """
        Create a mocpy object being the complemented of self

        Returns
        -------
        result : `~mocpy.AbstractMoc`
            the complemented moc

        """

        complement_interval = self._complement_interval()

        result = self.__class__.from_interval_set(complement_interval)
        return result

    def _complement_interval(self):
        from .interval_set import IntervalSet
        res = IntervalSet()

        interval_set = sorted(self._interval_set.intervals)

        if interval_set[0][0] > 0:
            res.add((0, interval_set[0][0]))

        last = interval_set[0][1]

        for itv in interval_set[1:]:
            res.add((last, itv[0]))
            last = itv[1]

        max_pix_order = self._get_max_pix()

        if last < max_pix_order:
            res.add((last, max_pix_order))

        return res

    def _get_max_pix(self):
        pass

    @staticmethod
    def _neighbour_pixels(hp, pix_arr):
        """
        Get all the pixels neighbours of ``pix_arr``

        Parameters
        ----------
        hp : `~astropy_healpix.HEALPix`
            the HEALPix context
        pix_arr : `~numpy.ndarray`
            the input array of pixels
        Returns
        -------
        neighbour_pix_arr : `~numpy.ndarray`
            an array of pixels containing the neighbours of the pixels in ``pix_arr``

        """
        neighbour_pix_arr = np.unique(hp.neighbours(pix_arr).ravel())
        # Remove negative pixel values returned by `~astropy_healpix.HEALPix.neighbours`
        neighbour_pix_arr = neighbour_pix_arr[np.where(neighbour_pix_arr >= 0)]
        return neighbour_pix_arr

    @classmethod
    def from_interval_set(cls, interval_set):
        """
        Create a MOC/TimeMoc from an `~mocpy.interval_set.IntervalSet` object (all intervals at deepest norder)

        Parameters
        ----------
        interval_set : `~mocpy.interval_set.IntervalSet`
            a set of intervals representing a MOC/TimeMoc

        Returns
        -------
        moc : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            the MOC/TimeMoc object reflecting ``interval_set``.

        """

        moc = cls()
        moc._interval_set = interval_set

        return moc

    @classmethod
    def from_json(cls, json_moc):
        """
        Create a MOC/TimeMoc from a dictionary of order : [pix] (.i.e json format)

        Parameters
        ----------
        json_moc : {str : [int]}
            The moc expressed as a dictionary of pixel lists indexed by their order

        Returns
        -------
        moc : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            the MOC/TimeMoc object reflecting ``json_moc``.

        """

        moc = cls()
        for order, i_pix_l in json_moc.items():
            moc.add_pix_list(order=int(order), i_pix_l=i_pix_l, nest=True)

        return moc

    @classmethod
    def from_uniq_interval_set(cls, uniq_is):
        """
        Create a MOC/TimeMoc from an IntervalSet of NUNIQ HEALPix pixels

        Parameters
        ----------
        uniq_is : `~mocpy.interval_set.IntervalSet`
            interval set of NUNIQ HEALPix pixels

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            the MOC/TimeMoc object reflecting ``uniq_is`` NUNIQ HEALPix interval set.

        """

        r = IntervalSet()
        rtmp = IntervalSet()
        last_order = 0

        intervals = uniq_is.intervals
        diff_order = AbstractMoc.HPY_MAX_NORDER
        shift_order = 2 * diff_order
        for interval in intervals:
            for j in range(interval[0], interval[1]):
                order, i_pix = utils.uniq2orderipix(j)

                if order != last_order:
                    r = r.union(rtmp)
                    rtmp.clear()
                    last_order = order
                    diff_order = AbstractMoc.HPY_MAX_NORDER - order
                    shift_order = 2 * diff_order

                rtmp.add((i_pix << shift_order, (i_pix + 1) << shift_order))

        r = r.union(rtmp)
        moc = cls.from_interval_set(r)
        return moc

    def uniq_pixels_iterator(self):
        """
        Generator giving the NUNIQ HEALPix pixels of the Moc/TimeMoc

        Returns
        -------
        uniq :
            the NUNIQ HEALPix pixels iterator

        """

        for uniq_iv in self.to_uniq_interval_set().intervals:
            for uniq in range(uniq_iv[0], uniq_iv[1]):
                yield uniq

    def to_uniq_interval_set(self):
        """
        Convert the default interval set form of the MOC/TimeMOC into an interval set of NUNIQ HEALPix pixels.

        Returns
        -------
        res : `~mocpy.interval_set.IntervalSet`
            the interval set of NUNIQ HEALPix pixels reflecting the current state of the MOC/TimeMoc.

        """

        r2 = IntervalSet(self._interval_set)
        r3 = IntervalSet()
        res = IntervalSet()

        intervals = self._interval_set.intervals
        if len(intervals) == 0:
            return res

        max_res_order = AbstractMoc.HPY_MAX_NORDER
        order = 0
        while not r2.empty():
            shift = int(2 * (max_res_order - order))
            ofs = (int(1) << shift) - 1
            ofs2 = int(1) << (2 * order + 2)

            r3.clear()
            for iv in r2.intervals:
                a = (int(iv[0]) + ofs) >> shift
                b = int(iv[1]) >> shift
                r3.add((a << shift, b << shift))
                res.add((a + ofs2, b + ofs2))

            if not r3.empty():
                r2 = r2.difference(r3)

            order += 1

        return res

    @classmethod
    def from_moc_fits_file(cls, path):
        """
        Load a MOC from a MOC fits file (i.e. a fits file in which pix are stored as a list of NUNIQ HEALPix pixels
        in a binary table HDU).

        It corresponds to the type of file in which the MOCs/TMOCs are stored when writing them
        in a fits file with the method `mocpy.AbstractMoc.write`.

        Parameters
        ----------
        path : str
            the path to the moc fits file

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            the mocpy object having as interval set the one stored in the fits file located at ``path``

        """

        with fits.open(path) as hdulist:
            # Case of a moc written in a fits file. Moc fits file stores the nuniq items in a
            # BinTableHDU usually at index 1
            if isinstance(hdulist[1], fits.hdu.table.BinTableHDU):
                interval_set = IntervalSet()
                # accessing directly recarray dramatically speed up the reading
                data = hdulist[1].data.view(np.recarray)
                for x in range(0, len(hdulist[1].data)):
                    interval_set.add(data[x][0])

                result = cls.from_uniq_interval_set(interval_set)
                return result

    def best_res_pixels_iterator(self):
        """
        Generator giving the pixels at the max order

        Returns
        -------
        val :
            the pixel iterator

        """

        factor = 4 ** (AbstractMoc.HPY_MAX_NORDER - self.max_order)
        for iv in self._interval_set.intervals:
            for val in range(iv[0] // factor, iv[1] // factor):
                yield val

    def _get_max_order_pix(self, keep_inside):
        """
        Get a boolean numpy array of size the number of pix in the max_order of the MOC object
        The ith element of this array equals to True if the corresponding pix is in the MOC for this order
        If it is not, the element is set to False

        Parameters
        ----------
        keep_inside : bool
            ``keep_inside`` boolean value is associated with pix that are in the MOC

        Returns
        -------
        result : `~numpy.ndarray`
            boolean array telling which pix at the max_order are in the MOC

        """

        from astropy_healpix.healpy import nside2npix
        max_order = self.max_order
        n_side = 2 ** max_order

        result = np.zeros(nside2npix(n_side), dtype=bool)
        for val in self.best_res_pixels_iterator():
            result[val] = True

        if not keep_inside:
            result = np.logical_not(result)

        return result

    def add_fits_header(self, tbhdu):
        """
        This method must be implemented in each class derived from AbstractMoc

        tbhdu : `~astropy.fits.BinTableHDU`
            fits HDU binary table

        """
        pass

    @staticmethod
    def _to_json(uniq_arr):
        """
        Serialize a mocpy object (array of uniq) to json

        Parameters
        ----------
        uniq_arr : `~numpy.ndarray`
            the array of uniq reprensenting the mocpy object to serialize

        Returns
        -------
        result_json : {str : [int]}
            a dictionary of pixel list each indexed by their order

        """

        result_json = {}

        order_arr, ipix_arr = utils.uniq2orderipix(uniq_arr)
        min_order = order_arr[0]
        max_order = order_arr[-1]

        for order in range(min_order, max_order+1):
            pix_index = np.where(order_arr == order)[0]
            if pix_index.size:
                # there are pixels belonging to the current order
                ipix_order_arr = ipix_arr[pix_index]
                result_json[str(order)] = ipix_order_arr.tolist()

        return result_json

    def _to_fits(self, uniq_arr, optional_kw_dict=None):
        """
        Serialize a mocpy object (array of uniq) to a fits format

        Parameters
        ----------
        uniq_arr : `~numpy.ndarray`
            the array of uniq representing the mocpy object to serialize
        optional_kw_dict : dict
            optional keywords arguments added to the fits header

        Returns
        -------
        thdulist : `~astropy.io.fits.HDUList`
            the fits serialization of the MOC/TimeMoc object

        """

        moc_order = self.max_order
        if moc_order <= 13:
            fits_format = '1J'
        else:
            fits_format = '1K'

        tbhdu = fits.BinTableHDU.from_columns(
            fits.ColDefs([fits.Column(name='UNIQ', format=fits_format, array=uniq_arr)]))
        tbhdu.header['PIXTYPE'] = 'HEALPIX'
        tbhdu.header['ORDERING'] = 'NUNIQ'
        self.add_fits_header(tbhdu)
        tbhdu.header['MOCORDER'] = moc_order
        tbhdu.header['MOCTOOL'] = 'MOCPy'
        if optional_kw_dict:
            for key in optional_kw_dict:
                tbhdu.header[key] = optional_kw_dict[key]

        thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
        return thdulist

    def write(self, path=None, format='fits', optional_kw_dict=None, write_to_file=False):
        """
        Serialize a MOC/TimeMoc object.

        Possibility to write it to a file at ``path``. Format can be 'fits' or 'json',
        though only the fits format is officially supported by the IVOA.

        Parameters
        ----------
        path : str, optional
            path to save the MOC object. The mocpy is written to path only if ``serialize`` is False. None by default
        format : str, optional
            format in which the mocpy object will be serialized. Constraint to takes its value
            among "fits" or "json". By default, ``format`` is set to "fits".
        optional_kw_dict : {str, _}, optional
            optional dictionary keywords for the header of the fits file. Only used if ``format`` is "fits"
        write_to_file : bool, optional
            Set to False by default. In this case, this method does not write to a file but returns the serialized form
            of the MOC/TimeMoc object to the user. If you want to write to a file

        Returns
        -------
        result : a `~astropy.io.fits.HDUList` if ``format`` is set to "fits" or {str, [int]} otherwise
            The serialization of the MOC/TimeMoc object

        """

        formats = ('fits', 'json')
        if format not in formats:
            raise ValueError('format should be one of %s' % (str(formats)))

        self._ensure_consistency()
        uniq_l = []

        for uniq in self.uniq_pixels_iterator():
            uniq_l.append(uniq)

        uniq_arr = np.array(uniq_l)

        if format == 'fits':
            result = self._to_fits(uniq_arr=uniq_arr,
                                   optional_kw_dict=optional_kw_dict)
            if write_to_file:
                result.writeto(path, overwrite=True)
        else:
            # json format serialization
            result = self.__class__._to_json(uniq_arr=uniq_arr)
            if write_to_file:
                import json
                with open(path, 'w') as h:
                    h.write(json.dumps(result, sort_keys=True, indent=2))

        return result

    def degrade_to_order(self, new_order):
        """
        Degrade self to a mocpy object at max_order being equal to ``new_order``

        Parameters
        ----------
        new_order : int

        Returns
        -------
        moc : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            the res decreased mocpy object

        """

        shift = 2 * (AbstractMoc.HPY_MAX_NORDER - new_order)
        ofs = (int(1) << shift) - 1
        mask = ~ofs
        adda = int(0)
        addb = ofs
        iv_set = IntervalSet()
        for iv in self._interval_set.intervals:
            a = (iv[0] + adda) & mask
            b = (iv[1] + addb) & mask
            if b > a:
                iv_set.add((a, b))

        moc = self.__class__.from_interval_set(iv_set)
        moc._order = new_order
        return moc
