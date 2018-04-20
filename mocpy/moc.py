#!/usr/bin/env python
# -*- coding: utf-8 -*

"""

moc.py:
  functions to read/write and manipulate MOCs

"""

__author__ = "Thomas Boch"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

try:
    set
except NameError:
    from sets import Set as set

import sys

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS
from astropy.io import fits
from astropy import wcs

from astropy_healpix import HEALPix
from astropy_healpix.healpy import ang2pix
from astropy_healpix.healpy import nside2npix

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from .abstract_moc import AbstractMoc
from .interval_set import IntervalSet
from . import utils

if sys.version > '3':
    long = int

# Python 3 support
try:
    xrange
except NameError:
    xrange = range


class MOC(AbstractMoc):
    VIZ_TABLE_MOC_ROOT_URL = ''
    VIZ_CAT_MOC_ROOT_URL = ''
    
    def __init__(self):
        AbstractMoc.__init__(self)

    def add_position(self, ra, dec, max_norder):
        """
        add the HEALPix bin containing the (ra, dec) position
        """
        theta, phi = utils.radec2thetaphi(ra, dec)
        i_pix = ang2pix(2**max_norder, theta, phi, nest=True)

        if isinstance(i_pix, np.ndarray):
            self.add_pix_list(order=max_norder, i_pix_l=i_pix)
        else:
            self.add_pix(max_norder, i_pix)

    def filter_table(self, table, ra_column, dec_column, keep_inside=True):
        return self._filter(table,
                            keep_inside,
                            None,
                            ra_column,
                            dec_column)

    def _get_pix(self, row_values_l, n_side, format=None):
        assert len(row_values_l) == 2, ValueError('Cannot filter following non spatial columns type.'
                                                  'ra and dec column are required')
        ra = row_values_l[0]
        dec = row_values_l[1]

        theta, phi = utils.radec2thetaphi(ra, dec)
        return ang2pix(n_side, theta, phi, nest=True)

    @classmethod
    def __from_fits_image(cls, path, moc_order):
        assert moc_order, ValueError('An order must be specified when'
                                     ' building a moc from a fits image.')
        # load the image data
        header = fits.getheader(path)
        height = header['NAXIS2']
        width = header['NAXIS1']

        # use wcs from astropy to locate the image in the world coordinates
        w = wcs.WCS(header)

        d_pix = 1
        x, y = np.mgrid[-0.25:(width + 1.25):d_pix, -0.25:(height + 1.25):d_pix]

        pix_crd = np.dstack((x.ravel(), y.ravel()))[0]

        world_pix_crd = w.wcs_pix2world(pix_crd, 1)
        hp = HEALPix(nside=(1 << moc_order), order='nested', frame=ICRS())

        i_pix_l = hp.lonlat_to_healpix(world_pix_crd[:, 0] * u.deg,
                                       world_pix_crd[:, 1] * u.deg)
        # remove doubles
        i_pix_l = list(set(i_pix_l))
        moc = MOC()
        moc.add_pix_list(order=moc_order, i_pix_l=i_pix_l)
        # this will be consistent when one will do operations on the moc (union, inter, ...) or
        # simply write it to a fits or json file

        return moc

    @classmethod
    def from_file(cls, path, moc_order=None):
        """
        Load a moc from a fits file (image or binary table)

        :param path: the path to the fits file
        :param moc_order: the max order in case one wants to create a moc from a fits image
        :return: a moc object corresponding to the passed fits file
        """
        with fits.open(path) as hdulist:
            if isinstance(hdulist[1], fits.hdu.table.BinTableHDU):
                # AbstractMoc returns a IntervalSet made of uniq
                return MOC.from_uniq_interval_set(AbstractMoc.from_file(hdulist=hdulist))

            assert isinstance(hdulist[1], fits.hdu.ImageHDU), ValueError('Cannot extract the moc'
                                                                         ' from the {0:s} fits file'.format(path))
            return MOC.__from_fits_image(path=path, moc_order=moc_order)

    @classmethod
    def from_fits_images(cls, path_l, moc_order):
        """
        Load a moc from a set of fits images

        :param path_l: the list of path where the fits image are located
        :param moc_order: the order one wants to build the moc
        :return: the union of all the moc from path_l
        """
        moc = MOC()
        for path in path_l:
            current_moc = MOC.from_file(path=path, moc_order=moc_order)
            moc = moc.union(current_moc)

        return moc

    @classmethod
    def from_vizier_table(cls, table_id, nside=256):
        """
        return the MOC for a given VizieR table
        """
        if nside not in (8, 16, 32, 64, 128, 256, 512):
            raise Exception('Bad value for nside')

        return cls.from_ivorn('ivo://CDS/' + table_id, nside)

    MOC_SERVER_ROOT_URL = 'http://alasky.unistra.fr/MocServer/query'

    @classmethod
    def from_ivorn(cls, ivorn, nside=256):
        """
        return the MOC for a given IVORN
        """
        try:
            from urllib import urlencode
        except ImportError:
            from urllib.parse import urlencode

        return cls.from_url('%s?%s' % (MOC.MOC_SERVER_ROOT_URL,
                                       urlencode({
                                           'ivorn': ivorn,
                                           'get': 'moc',
                                           'order': int(np.log2(nside))
                                       })))

    @classmethod
    def from_url(cls, url):
        from astropy.utils.data import download_file

        path = download_file(url, show_progress=False, timeout=60)
        return cls.from_file(path)

    @classmethod
    def from_table(cls, table, ra_column, dec_column, moc_order):
        """
        Create a MOC from a astropy.table.Table
        The user has to specify the columns holding ra and dec (in ICRS)
        """
        if not isinstance(ra_column, str) or not isinstance(dec_column, str):
            raise TypeError

        moc = MOC()
        moc._order = moc_order

        moc.add_position(table[ra_column], table[dec_column], moc_order)

        return moc

    @classmethod
    def from_coo_list(cls, skycoords, max_norder):
        """
        Create a MOC from a list of SkyCoord
        """
        if not isinstance(skycoords, SkyCoord):
            raise TypeError

        moc = MOC()
        moc._order = max_norder

        moc.add_position(skycoords.icrs.ra.deg, skycoords.icrs.dec.deg, max_norder)

        return moc

    @property
    def sky_fraction(self):
        """
        return the sky fraction (between 0 and 1) covered by the MOC
        """
        nb_pix_filled = 0
        for val in self.best_res_pixels_iterator():
            nb_pix_filled += 1

        return nb_pix_filled / float(3 << (2*(self.max_order + 1)))

    def query_simbad(self, max_rows=10000):
        """
        query a view of SIMBAD data
        for SIMBAD objects in the coverage of the MOC instance
        """
        return self._query('SIMBAD', max_rows)

    def query_vizier_table(self, table_id, max_rows=10000):
        """
        query a VizieR table
        for sources in the coverage of the MOC instance
        """
        return self._query(table_id, max_rows)

    def _query(self, resource_id, max_rows):
        """
        internal method to query Simbad or a VizieR table
        for sources in the coverage of the MOC instance
        """
        import requests
        import tempfile
        import os

        if max_rows is not None and max_rows >= 0:
            max_rows_str = str(max_rows)
        else:
            max_rows_str = str(9999999999)

        tmp_moc = tempfile.NamedTemporaryFile(delete=False)
        self.write(tmp_moc.name)
        r = requests.post('http://cdsxmatch.u-strasbg.fr/QueryCat/QueryCat',
                          data={'mode': 'mocfile',
                                'catName': resource_id,
                                'format': 'votable',
                                'limit': max_rows_str},
                          files={'moc': open(tmp_moc.name, 'rb')},
                          headers={'User-Agent': 'MOCPy'},
                          stream=True)
        
        tmp_vot = tempfile.NamedTemporaryFile(delete=False)
        with open(tmp_vot.name, 'w') as h:
            for line in r.iter_lines():
                if line:
                    h.write(line.decode(r.encoding)+'\n')

        from astropy.io.votable import parse_single_table
        table = parse_single_table(tmp_vot.name).to_table()

        # finally delete temp files
        os.unlink(tmp_moc.name)
        os.unlink(tmp_vot.name)

        return table

    """
    Plot a spatial moc using matplotlib and healpy 
    """
    # TODO : remove all healpy dependencies (e.g. mollview and graticule)
    def degrade_to_order(self, new_order):
        shift = 2 * (MOC.HPY_MAX_NORDER - new_order)
        ofs = (long(1) << shift) - 1
        mask = ~ofs
        adda = long(0)
        addb = ofs
        iv_set = IntervalSet()
        for iv in self._interval_set.intervals:
            a = (iv[0] + adda) & mask
            b = (iv[1] + addb) & mask
            if b > a:
                iv_set.add((a, b))

        m = MOC.from_interval_set(iv_set)
        m._order = new_order
        return m

    def plot(self, title='MOC', coord='C'):
        plot_order = 8
        if self.max_order > plot_order:
            plotted_moc = self.degrade_to_order(plot_order)
        else:
            plotted_moc = self

        num_pixels_map = 768
        delta = 2 * np.pi / num_pixels_map

        x = np.arange(-np.pi, np.pi, delta)
        y = np.arange(-np.pi/2, np.pi/2, delta)
        lon_rad, lat_rad = np.meshgrid(x, y)

        m = np.zeros(nside2npix(2 ** plotted_moc.max_order))
        for val in plotted_moc.best_res_pixels_iterator():
            m[val] = 1

        hp = HEALPix(nside=(1 << plotted_moc.max_order), order='nested', frame=ICRS())

        pix_map = hp.lonlat_to_healpix(lon_rad * u.rad, lat_rad * u.rad)
        map = np.flip(m[pix_map], axis=1)

        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111, projection='mollweide')
        from matplotlib.colors import LinearSegmentedColormap
        color_map = LinearSegmentedColormap.from_list('w2r', ['#eeeeee', '#aa0000'])
        color_map.set_under('w')
        color_map.set_bad('gray')

        ax.pcolormesh(x, y, map, cmap=color_map, vmin=0, vmax=1)
        ax.tick_params(labelsize=14, labelcolor='#000000')
        plt.title('moc plot:')
        plt.grid(True, linestyle='--', linewidth=1, color='#555555')

        plt.show()
