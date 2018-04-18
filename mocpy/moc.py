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
from astropy.table import Table
from astropy.coordinates import SkyCoord, ICRS
from astropy.io import fits
from astropy import wcs

from astropy_healpix import HEALPix
from astropy_healpix.healpy import ang2pix


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
        ipix = ang2pix(2**max_norder, theta, phi, nest=True)

        if isinstance(ipix, np.ndarray):
            self.add_pix_list(order=max_norder, i_pix_l=ipix)
        else:
            self.add_pix(max_norder, ipix)

    @classmethod
    def __from_fits_image(cls, path, max_order=None):
        assert max_order, ValueError('An order must be specified when'
                                     ' building a moc from a fits image.')
        # load the image data
        data_image = fits.getdata(path, ext=0)
        header = fits.getheader(path)
        height = header['NAXIS2']
        width = header['NAXIS1']
        print(height, width)

        # use wcs from astropy to locate the image in the world coordinates
        w = wcs.WCS(header)

        d_pix = 1
        X, Y = np.mgrid[0:(width + 1):d_pix, 0:(height + 1):d_pix]

        pix_crd = np.dstack((X.ravel(), Y.ravel()))[0]

        world_pix_crd = w.wcs_pix2world(pix_crd, 1)
        hp = HEALPix(nside=(1 << max_order), order='nested', frame=ICRS())

        print(world_pix_crd)

        ipix_l = hp.lonlat_to_healpix(world_pix_crd[:, 0] * u.deg,
                                      world_pix_crd[:, 1] * u.deg)
        # remove doubles
        ipix_l = list(set(ipix_l))
        moc = MOC()
        moc.add_pix_list(order=max_order, i_pix_l=ipix_l)

        return moc

    @classmethod
    def from_file(cls, path, max_order=None):
        """
        Load a moc from a fits file

        :param path: the path to the fits file
        :param max_order: the max order in case one wants to create a moc from a fits image
        :return: a moc object corresponding to the passed fits file
        """
        with fits.open(path) as hdulist:
            if isinstance(hdulist[1], fits.hdu.table.BinTableHDU):
                return AbstractMoc.from_file(hdulist=hdulist)

            assert isinstance(hdulist[1], fits.hdu.ImageHDU), ValueError('Cannot extract the moc'
                                                                         ' from the {0:s} fits file'.format(path))
            return MOC.__from_fits_image(path=path, max_order=max_order)

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

        return nb_pix_filled / float(12 * 4**self.max_order)

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

        if max_rows is not None and max_rows>=0:
            max_rows_str = str(max_rows)
        else:
            max_rows_str = str(9999999999)

        tmp_moc = tempfile.NamedTemporaryFile(delete = False)
        self.write(tmp_moc.name)
        r = requests.post('http://cdsxmatch.u-strasbg.fr/QueryCat/QueryCat', data={'mode': 'mocfile' , 'catName': resource_id, 'format': 'votable', 'limit': max_rows_str}, files={'moc': open(tmp_moc.name, 'rb')}, headers={'User-Agent': 'MOCPy'}, stream=True)
        
        tmp_vot = tempfile.NamedTemporaryFile(delete = False)
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

    def best_res_pixels_iterator(self):
        factor = 4**(MOC.HPY_MAX_NORDER - self.max_order)
        for iv in self._interval_set.intervals:
            for val in xrange(iv[0] // factor, iv[1] // factor):
                yield val
    
    def filter_table(self, table, ra_column, dec_column, keep_inside=True):
        """
        Filter an astropy.table.Table to keep only rows inside (or outside) the MOC instance
        Return the (newly created) filtered Table
        """
        filtered_table = Table()
        
        kept_rows = []
        pixels_best_res = set()
        for val in self.best_res_pixels_iterator():
            pixels_best_res.add(val)

        max_order = self.max_order
        nside = 2**max_order
        for row in table:
            theta, phi = utils.radec2thetaphi(row[ra_column], row[dec_column])
            ipix = ang2pix(nside, theta, phi, nest=True)
            if (ipix in pixels_best_res) == keep_inside:
                kept_rows.append(row)

        if len(kept_rows) == 0:
            return Table(names=table.colnames)
        else:
            return Table(rows=kept_rows, names=table.colnames)


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
        import matplotlib.pyplot as plt

        import numpy as np
        from astropy_healpix import HEALPix
        from astropy.coordinates import ICRS, Galactic

        if self.max_order > 4:
            plotted_moc = self.degrade_to_order(4)
        else:
            plotted_moc = self

        hp = HEALPix(nside=(1 << plotted_moc.max_order), order='nested', frame=Galactic())
        pix_best_res_l = list(plotted_moc.best_res_pixels_iterator())

        lon, lat = hp.boundaries_lonlat(pix_best_res_l, step=1)
        print(lon.to(u.deg).value)
        print(lat.to(u.deg).value)
        deg_to_rad = np.pi / 180
        lon = (lon.to(u.deg).value - 180) * deg_to_rad
        lat = lat.to(u.deg).value * deg_to_rad

        print('max', lon.max())

        plt.figure()
        ax = plt.subplot(111, projection="mollweide")
        for i in range(len(lon)):
            ax.fill(lon[i], lat[i], 'r')

        plt.title("Mollweide")
        plt.grid(True)

        plt.show()

    '''
    def plot(self, title='MOC', coord='C'):
        """
        plot current instance using matplotlib
        coord can be 'C', 'G' or 'E' (respectively celestial, galactic or ecliptic)
        """
        from matplotlib import pyplot as plt
        import healpy as hp
        # degrade to NORDER 8 if norder is greater
        if self.max_order > 8:
            plotted_moc = self.degrade_to_order(8)
        else:
            plotted_moc = self
        import pdb; pdb.set_trace()
        m = np.zeros(hp.nside2npix(2 ** plotted_moc.max_order))
        for val in plotted_moc.best_res_pixels_iterator():
            m[val] = 1

        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list('w2r', ['#ffffff', '#ff0000'])
        cmap.set_under('w') 
        cmap.set_bad('gray')
        hp.mollview(m, nest=True, coord=['C', coord], title=title, cbar=False, cmap=cmap)
        hp.graticule()
        plt.show()
    '''