# -*- coding: utf-8 -*
from __future__ import absolute_import, division, print_function, unicode_literals
from .py23_compat import range, int

import numpy as np

from astropy.io import fits
from astropy.table import Table

from .interval_set import IntervalSet
from . import utils

__author__ = "Thomas Boch, Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, matthieu.baumann@astro.unistra.fr"


class AbstractMOC:
    """
    Basic functions for manipulating MOCs.
    """
    HPY_MAX_NORDER = 29

    def __init__(self, interval_set=None):
        interval = IntervalSet() if interval_set is None else interval_set
        self._interval_set = interval
        # Must be overridden in subclasses
        self._fits_header_keywords = None

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
        if not isinstance(another_moc, AbstractMOC):
            raise TypeError('Cannot compare an AbstractMOC with a {0}'.format(type(another_moc)))

        return self._interval_set == another_moc._interval_set

    @property
    def max_order(self):
        """
        Returns the max depth of the MOC.

        The max depth of a MOC is the depth of the smallest HEALPix cells that are contained
        in the MOC.
        """
        # TODO: cache value
        combo = int(0)
        for iv in self._interval_set._intervals:
            combo |= iv[0] | iv[1]

        ret = AbstractMOC.HPY_MAX_NORDER - (utils.number_trailing_zeros(combo) // 2)
        if ret < 0:
            ret = 0

        return ret

    def intersection(self, another_moc, *args):
        """
        Intersection between self and other MOCs.

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMoc used for performing the intersection with self
        args : `~mocpy.abstract_moc.AbstractMOC`
            other MOCs

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self & ``moc``
        """
        interval_set = self._interval_set.intersection(another_moc._interval_set)
        for moc in args:
            interval_set = interval_set.intersection(moc._interval_set)

        return self.__class__(interval_set)

    def union(self, another_moc, *args):
        """
        Union between self and other MOCs.

        Parameters
        ----------
        another_moc : `mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMoc to bind to self
        args : `~mocpy.abstract_moc.AbstractMOC`
            other MOCs

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self | ``moc``
        """
        interval_set = self._interval_set.union(another_moc._interval_set)
        for moc in args:
            interval_set = interval_set.union(moc._interval_set)

        return self.__class__(interval_set)

    def difference(self, another_moc, *args):
        """
        Difference between self and other MOCs.

        Parameters
        ----------
        moc : `mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMoc to substract from self
        args : `~mocpy.abstract_moc.AbstractMOC`
            other MOCs

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self - ``moc``
        """
        interval_set = self._interval_set.difference(another_moc._interval_set)
        for moc in args:
            interval_set = interval_set.difference(moc._interval_set)

        return self.__class__(interval_set)

    def complement(self):
        """
        Return the complement of the MOC.

        Returns
        -------
        complement : `~mocpy.AbstractMoc`
            the complemented moc
        """
        complement_interval = self._complement_interval()
        return self.__class__(complement_interval)

    def _complement_interval(self):
        res = []
        intervals_l = sorted(self._interval_set._intervals.tolist())

        if intervals_l[0][0] > 0:
            res.append((0, intervals_l[0][0]))

        last = intervals_l[0][1]

        for itv in intervals_l[1:]:
            res.append((last, itv[0]))
            last = itv[1]

        max_pix_order = self._get_max_pix()

        if last < max_pix_order:
            res.append((last, max_pix_order))

        return IntervalSet.from_numpy_array(np.asarray(res))

    def _get_max_pix(self):
        pass

    @staticmethod
    def _neighbour_pixels(hp, ipix):
        """
        Get all the pixels neighbours of ``ipix``

        Parameters
        ----------
        hp : `~astropy_healpix.HEALPix`
            the HEALPix context
        ipix : `~numpy.ndarray`
            the input array of pixels
        Returns
        -------
        result : `~numpy.ndarray`
            an array of pixels containing the neighbours of the pixels in ``ipix``
        """
        neigh_ipix = np.unique(hp.neighbours(ipix).ravel())
        # Remove negative pixel values returned by `~astropy_healpix.HEALPix.neighbours`
        return neigh_ipix[np.where(neigh_ipix >= 0)]

    @classmethod
    def from_json(cls, json_moc):
        """
        Creates a MOC from a dictionary of HEALPix arrays indexed by their depth.

        Parameters
        ----------
        json_moc : {str : [int]}
            A dictionary of HEALPix arrays indexed by their order

        Returns
        -------
        moc : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            the MOC/TimeMoc object reflecting ``json_moc``.

        """
        intervals_arr = np.array([])
        for order, pix_l in json_moc.items():
            if len(pix_l) == 0:
                continue
            pix_arr = np.array(pix_l)
            p1 = pix_arr
            p2 = pix_arr + 1
            shift = 2 * (AbstractMOC.HPY_MAX_NORDER - int(order))

            itv_arr = np.vstack((p1 << shift, p2 << shift)).T
            if intervals_arr.size == 0:
                intervals_arr = itv_arr
            else:
                intervals_arr = np.vstack((intervals_arr, itv_arr))

        return cls(IntervalSet.from_numpy_array(intervals_arr))

    def _uniq_pixels_iterator(self):
        """
        Generator giving the NUNIQ HEALPix pixels of the Moc/TimeMoc

        Returns
        -------
        uniq :
            the NUNIQ HEALPix pixels iterator
        """
        intervals_uniq_l = IntervalSet.to_nuniq_interval_set(self._interval_set)._intervals
        for uniq_iv in intervals_uniq_l:
            for uniq in range(uniq_iv[0], uniq_iv[1]):
                yield uniq

    @classmethod
    def from_fits(cls, filename):
        """
        Load a MOC from a MOC fits file.

        It corresponds to the default type of file in which the MOCs/TMOCs are stored when using
        the method `mocpy.AbstractMoc.write`. 
        A fits file stores the list of NUNIQ HEALPix cells describing the MOC in a binary HDU table.

        Parameters
        ----------
        filename : str
            The path to the moc fits file

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            The resulting MOC.
        """
        table = Table.read(filename)

        intervals = np.vstack((table['UNIQ'], table['UNIQ']+1)).T

        nuniq_interval_set = IntervalSet.from_numpy_array(intervals)
        interval_set = IntervalSet.from_nuniq_interval_set(nuniq_interval_set)
        return cls(interval_set)

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
        min_order = np.min(order_arr[0])
        max_order = np.max(order_arr[-1])

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
        uniq_arr : `numpy.ndarray`
            the array of uniq representing the mocpy object to serialize
        optional_kw_dict : dict
            optional keywords arguments added to the fits header

        Returns
        -------
        thdulist : `astropy.io.fits.HDUList`
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
        tbhdu.header.update(self._fits_header_keywords)
        tbhdu.header['MOCORDER'] = moc_order
        tbhdu.header['MOCTOOL'] = 'MOCPy'
        if optional_kw_dict:
            for key in optional_kw_dict:
                tbhdu.header[key] = optional_kw_dict[key]

        thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
        return thdulist

    def serialize(self, format='fits', optional_kw_dict=None):
        """
        Serialize the MOC into a specific format.

        Parameters
        ----------
        format : str
            'fits' by default. Other choice possible is 'json'.
        optional_kw_dict : dict
            optional keywords arguments added to the fits header

        Returns
        -------
        result : `astropy.io.fits.HDUList` or json dictionary
            Serialized MOC.
        """
        formats = ('fits', 'json')
        if format not in formats:
            raise ValueError('format should be one of %s' % (str(formats)))

        uniq_l = []
        for uniq in self._uniq_pixels_iterator():
            uniq_l.append(uniq)

        uniq_arr = np.array(uniq_l)

        if format == 'fits':
            result = self._to_fits(uniq_arr=uniq_arr,
                                   optional_kw_dict=optional_kw_dict)
        else:
            # json format serialization
            result = self.__class__._to_json(uniq_arr=uniq_arr)

        return result

    def write(self, path, format='fits', optional_kw_dict=None):
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
        optional_kw_dict : optional
            optional dictionary keywords for the header of the fits file. Only used if ``format`` is "fits"
        write_to_file : bool, optional
            Set to False by default. In this case, this method does not write to a file but returns the serialized form
            of the MOC/TimeMoc object to the user. If you want to write to a file

        Returns
        -------
        result : `astropy.io.fits.HDUList` or JSON dict
            The serialization depending on the value of ``format``.
        """
        serialization = self.serialize(format=format, optional_kw_dict=optional_kw_dict)
        if format == 'fits':
            serialization.writeto(path, overwrite=True)
        else:
            import json
            with open(path, 'w') as h:
                h.write(json.dumps(serialization, sort_keys=True, indent=2))

    def degrade_to_order(self, new_order):
        """
        Degrade self to a mocpy object at max_order being equal to ``new_order``

        Parameters
        ----------
        new_order : int

        Returns
        -------
        moc : `mocpy.moc.MOC` or `mocpy.tmoc.TimeMoc`
            the res decreased mocpy object
        """
        shift = 2 * (AbstractMOC.HPY_MAX_NORDER - new_order)
        ofs = (int(1) << shift) - 1
        mask = ~ofs
        adda = int(0)
        addb = ofs
        iv_set = []

        for iv in self._interval_set._intervals:
            a = (iv[0] + adda) & mask
            b = (iv[1] + addb) & mask
            if b > a:
                iv_set.append((a, b))

        return self.__class__(IntervalSet.from_numpy_array(np.asarray(iv_set)))
