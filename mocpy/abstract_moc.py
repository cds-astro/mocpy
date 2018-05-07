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
        if not isinstance(another_moc, AbstractMoc):
            raise TypeError('Cannot compare an AbstractMoc with a {0}'.format(type(another_moc)))

        return self._interval_set.intervals == another_moc._interval_set.intervals

    def add_pix_list(self, order, i_pix_l, nest=True):
        """
        add a list of HEALPix pixel number at a given order
        to the current object. This uses numpy vectorization
        """
        add_pix_vect = np.vectorize(self.add_pix)
        add_pix_vect(order, i_pix_l, nest=nest)

    def add_pix(self, order, i_pix, nest=True):
        """
        add a given HEALPix pixel number to the current object

        """

        """
        # Not good for the vectorization
        self.__counter += 1
        # force consistency to prevent too large interval array
        if self.__counter == 1000:
            self._ensure_consistency()
            self.__counter = 0
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
        Force IntervalSet internal consistency
        """
        self._interval_set.intervals

    def intersection(self, another_moc):
        """
        intersection with another MOC
        """
        iv_set_intersection = self._interval_set.intersection(another_moc._interval_set)

        return self.__class__.from_interval_set(iv_set_intersection)

    def union(self, another_moc):
        """
        union with another MOC
        """
        iv_set_union = self._interval_set.union(another_moc._interval_set)

        return self.__class__.from_interval_set(iv_set_union)

    def difference(self, moc):
        """
        difference of self and moc
        :param moc: the moc to subtract from self
        :return: the result of the difference : self - moc
        """
        iv_set_difference = self._interval_set.difference(moc._interval_set)

        return self.__class__.from_interval_set(iv_set_difference)

    @classmethod
    def from_interval_set(cls, interval_set):
        """
        Create a MOC from an IntervalSet (all intervals at deepest norder)
        """
        moc = cls()
        moc._interval_set = interval_set

        return moc

    @classmethod
    def from_json(cls, json_moc):
        moc = cls()
        for order, i_pix_l in json_moc.items():
            moc.add_pix_list(order=int(order), i_pix_l=i_pix_l, nest=True)

        return moc

    @classmethod
    def from_uniq_interval_set(cls, uniq_is):
        """
        Create a MOC from an IntervalSet of NUNIQ HEALPix pixels

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
        return cls.from_interval_set(r)

    def uniq_pixels_iterator(self):
        for uniq_iv in self.to_uniq_interval_set().intervals:
            for uniq in range(uniq_iv[0], uniq_iv[1]):
                yield uniq

    def to_uniq_interval_set(self):
        """
        Return an IntervalSet in NUNIQ norder
        reflecting the current state of the MOC

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
        Load a MOC from a MOC fits file (i.e. a fits file in which pix are stored as a list of nuniq
        in a binary table HDU).

        It corresponds to the type of file in which the MOCs/TMOCs are stored when writing them
        in a fits file with the method `mocpy.AbstractMoc.write`.

        Parameters
        ----------
        path : str
            the path to the moc fits file

        Returns
        -------
            a mocpy.MOC/TimeMoc object

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
                return cls.from_uniq_interval_set(interval_set)

    def best_res_pixels_iterator(self):
        factor = 4 ** (AbstractMoc.HPY_MAX_NORDER - self.max_order)
        for iv in self._interval_set.intervals:
            for val in range(iv[0] // factor, iv[1] // factor):
                yield val

    def _get_max_order_pix(self, keep_inside):
        from astropy_healpix.healpy import nside2npix
        max_order = self.max_order
        n_side = 2 ** max_order

        m = np.zeros(nside2npix(n_side), dtype=bool)
        for val in self.best_res_pixels_iterator():
            m[val] = True

        if not keep_inside:
            m = np.logical_not(m)

        return m

    def add_fits_header(self, tbhdu):
        """ This method must be implemented in each class derived from AbstractMoc """
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

        m = self.__class__.from_interval_set(iv_set)
        m._order = new_order
        return m
