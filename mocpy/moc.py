#!/usr/bin/env python
# -*- coding: utf-8 -*

"""
moc.py

functions to read/write and manipulate MOCs

"""

from __future__ import absolute_import, division, print_function
from . import py23_compat

import requests
import tempfile
import os
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
from astropy import wcs
from astropy_healpix import HEALPix
from astropy_healpix.healpy import ang2pix
from astropy_healpix.healpy import nside2npix
import matplotlib.pyplot as plt

from .abstract_moc import AbstractMoc
from . import utils

__author__ = "Thomas Boch, Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, matthieu.baumann@astro.unistra.fr"


class MOC(AbstractMoc):
    VIZ_TABLE_MOC_ROOT_URL = ''
    VIZ_CAT_MOC_ROOT_URL = ''
    
    def __init__(self):
        AbstractMoc.__init__(self)

    def add_position(self, ra, dec, max_norder):
        """
        Add the HEALPix bin containing the (ra, dec) position

        Parameters
        ----------
        ra : `~astropy.coordinates.angles.Longitude`
            the longitude of the pixel to add
        dec : `~astropy.coordinates.angles.Latitude`
            the latitude of the pixel to add
        max_norder : int
            the moc order resolution

        """

        theta, phi = utils.radec2thetaphi(ra, dec)
        i_pix = ang2pix(2**max_norder, theta, phi, nest=True)

        if isinstance(i_pix, np.ndarray):
            self.add_pix_list(order=max_norder, i_pix_l=i_pix)
        else:
            self.add_pix(max_norder, i_pix)

    def filter_table(self, table, ra_column, dec_column, keep_inside=True):
        """
        Filter an `~astropy.table.Table` to keep only rows inside (or outside) the mocpy object instance

        Parameters
        ----------
        table : `~astropy.table.Table`
        ra_column : str
            name of the RA column to consider in the table
        dec_column : str
            name of the DEC column to consider in the table
        keep_inside : bool, optional
            True by default. In this case, the filtered table contains only observations that are located into
            the mocpy object. If ``keep_inside`` is False, the filtered table contains all observations lying outside
            the mocpy object.

        Returns
        -------
        table : `~astropy.table.Table`
            The (newly created) filtered Table

        """

        m = self._get_max_order_pix(keep_inside=keep_inside)
        theta, phi = utils.radec2thetaphi(table[ra_column], table[dec_column])
        pix_arr = ang2pix(2 ** self.max_order, theta, phi, nest=True)

        filtered_rows = m[pix_arr]
        return table[filtered_rows]

    def _get_max_pix(self):
        return 3*(2**60)

    def add_neighbours(self):
        """
        Add all the pixels at max order in the neighbourhood of the moc

        """
        from astropy_healpix import HEALPix

        hp = HEALPix(nside=(1 << self.max_order), order='nested')

        pix_arr = np.array(list(self.best_res_pixels_iterator()))
        neighbour_pix_arr = AbstractMoc._neighbour_pixels(hp, pix_arr)

        neighbour_pix_arr = np.setdiff1d(neighbour_pix_arr, pix_arr)

        factor = 4 ** (self.HPY_MAX_NORDER - self.max_order)
        for pix in neighbour_pix_arr:
            self._interval_set.add((pix * factor, (pix + 1) * factor))

    def remove_neighbours(self):
        """
        Remove all the pixels at max order located at the bound of the moc

        """
        from astropy_healpix import HEALPix

        hp = HEALPix(nside=(1 << self.max_order), order='nested')

        pix_arr = np.array(list(self.best_res_pixels_iterator()))

        neighbour_pix_arr = AbstractMoc._neighbour_pixels(hp, pix_arr)

        only_neighbour_arr = np.setxor1d(neighbour_pix_arr, pix_arr)

        bound_pix_arr = AbstractMoc._neighbour_pixels(hp, only_neighbour_arr)

        result_pix_arr = np.setdiff1d(pix_arr, bound_pix_arr)

        factor = 4 ** (self.HPY_MAX_NORDER - self.max_order)
        self._interval_set.clear()
        for pix in result_pix_arr:
            self._interval_set.add((pix * factor, (pix + 1) * factor))

    @classmethod
    def from_image(cls, header, hdu, moc_order, mask_arr=None):
        """
        Create a `~mocpy.moc.MOC` from an image stored as a fits file

        Parameters
        ----------
        header : `~astropy.io.fits.Header`
            fits header containing all the info of where the image is located (position, size, etc...)
        moc_order : int
            the moc resolution
        mask_arr : `~numpy.ndarray`, optional
            a 2D boolean array of the same size of the image where pixels having the value 1 are part of
            the survey and pixels having the value 0 are not part of the survey.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            the MOC object loaded from the ``mask_arr`` and ``header`` extracted from the image

        """
        # load the image data
        height = header['NAXIS2']
        width = header['NAXIS1']

        # use wcs from astropy to locate the image in the world coordinates
        w = wcs.WCS(header)

        try:
            if header['BLANK']:
                mask_arr = np.isfinite(hdu)
        except KeyError:
            pass

        if mask_arr is not None:
            # We have an array of pixels that are part of of survey
            y, x = np.where(mask_arr)
            pix_crd = np.dstack((x, y))[0]
        else:
            # If we do not have a mask array we create the moc of all the image
            #
            step_pix = 1
            """
            Coords returned by wcs_pix2world method correspond to pixel centers. We want to retrieve the moc pix
            crossing the borders of the image so we have to add 1/2 to the pixels coords before computing the lonlat.
            
            The step between two pix_crd is set to `step_pix` but can be diminished to have a better precision at the 
            borders so that all the image is covered (a too big step does not retrieve all
            the moc pix crossing the borders of the image).
            """
            x, y = np.mgrid[0.5:(width + 0.5 + step_pix):step_pix, 0.5:(height + 0.5 + step_pix):step_pix]
            pix_crd = np.dstack((x.ravel(), y.ravel()))[0]

        world_pix_crd = w.wcs_pix2world(pix_crd, 1)

        hp = HEALPix(nside=(1 << moc_order), order='nested', frame=ICRS())
        i_pix_l = hp.lonlat_to_healpix(lon=world_pix_crd[:, 0] * u.deg,
                                       lat=world_pix_crd[:, 1] * u.deg)
        # remove doubles
        i_pix_l = np.unique(i_pix_l)
        moc = MOC()
        moc.add_pix_list(order=moc_order, i_pix_l=i_pix_l)
        # this will be consistent when one will do operations on the moc (union, inter, ...) or
        # simply write it to a fits or json file

        return moc

    @classmethod
    def from_fits_images(cls, path_l, moc_order):
        """
        Load a moc from a set of fits images

        Parameters
        ----------
        path_l : [str]
            the path list where the fits image are located
        moc_order : int
            moc resolution

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            the union of all the moc from path_l

        """

        from astropy.io import fits
        moc = MOC()
        for path in path_l:
            header = fits.getheader(path)
            current_moc = MOC.from_image(header=header, moc_order=moc_order)
            moc = moc.union(current_moc)

        return moc

    @classmethod
    def from_vizier_table(cls, table_id, nside=256):
        """
        Create a `~mocpy.moc.MOC` object from a VizieR table

        Parameters
        ----------
        table_id : int
            table index
        nside : int, optional
            256 by default

        Returns
        -------
        result : `~mocpy.moc.MOC`
            the created moc

        """

        nside_possible_values = (8, 16, 32, 64, 128, 256, 512)
        if nside not in nside_possible_values:
            raise ValueError('Bad value for nside. Must be in {0}'.format(nside_possible_values))

        result = cls.from_ivorn('ivo://CDS/' + table_id, nside)
        return result

    MOC_SERVER_ROOT_URL = 'http://alasky.unistra.fr/MocServer/query'

    @classmethod
    def from_ivorn(cls, ivorn, nside=256):
        """
        Create a `~mocpy.moc.MOC` object from a given ivorn

        Parameters
        ----------
        ivorn : str
        nside : int, optional
            256 by default

        Returns
        -------
        result : `~mocpy.moc.MOC`
            the created moc

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
        """
        Create a `~mocpy.moc.MOC` object from a given url

        Parameters
        ----------
        url : str
            the url where the fits file representing a moc is available

        Returns
        -------
        result : `~mocpy.moc.MOC`
            the created moc

        """

        from astropy.utils.data import download_file

        path = download_file(url, show_progress=False, timeout=60)
        result = cls.from_moc_fits_file(path)
        return result

    @classmethod
    def from_table(cls, table, ra_column, dec_column, moc_order):
        """
        Create a `~mocpy.moc.MOC` object from a `~astropy.table.Table`

        The user has to specify the columns holding the ra and dec (in ICRS) data in ``table``

        Parameters
        ----------
        table : `~astropy.table.Table`
            the observations astropy table
        ra_column : str
            the name of the column referring to the right ascension of the observations
        dec_column : str
            the name of the column referring to the declination of the observations
        moc_order : int
            the max order of the MOC that will be created

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            the created moc

        """

        if not isinstance(ra_column, str) or not isinstance(dec_column, str):
            raise TypeError('`ra_column` whose type is {0} or `dec_column` whose type is {1} must be of string type'
                            .format())

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

        nb_pix_filled = len(list(self.best_res_pixels_iterator()))
        return nb_pix_filled / float(3 << (2*(self.max_order + 1)))

    def query_simbad(self, max_rows=10000):
        """
        query a view of SIMBAD data for SIMBAD objects in the coverage of the MOC instance

        """

        return self._query('SIMBAD', max_rows)

    def query_vizier_table(self, table_id, max_rows=10000):
        """
        query a VizieR table for sources in the coverage of the MOC instance

        """

        return self._query(table_id, max_rows)

    def _query(self, resource_id, max_rows):
        """
        internal method to query Simbad or a VizieR table
        for sources in the coverage of the MOC instance

        """

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

        from six import BytesIO
        tmp_vot = BytesIO()
        tmp_vot.write(r.content)

        from astropy.io.votable import parse_single_table
        table = parse_single_table(tmp_vot).to_table()

        # finally delete temp files
        os.unlink(tmp_moc.name)

        return table

    def add_fits_header(self, tbhdu):
        tbhdu.header['COORDSYS'] = ('C', 'reference frame (C=ICRS)')

    def plot(self, title='MOC', coord='C'):
        """
        Plot the MOC object in a mollweide view

        This method uses matplotlib.

        Parameters
        ----------
        title : str
            the title of the plot
        coord : str
            type of coord (ICRS, Galactic, ...) in which the moc pix will be plotted.
            only ICRS coordinates are supported for the moment.
            #TODO handle Galactic coordinates

        """

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
        z = np.flip(m[pix_map], axis=1)

        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111, projection="mollweide")
        ax.set_xticklabels(['150°', '120°', '90°', '60°', '30°', '0°', '330°', '300°', '270°', '240°', '210°', '180°'])

        from matplotlib.colors import LinearSegmentedColormap
        color_map = LinearSegmentedColormap.from_list('w2r', ['#eeeeee', '#aa0000'])
        color_map.set_under('w')
        color_map.set_bad('gray')

        ax.pcolormesh(x, y, z, cmap=color_map, vmin=0, vmax=1)
        ax.tick_params(labelsize=14, labelcolor='#000000')
        plt.title(title)
        plt.grid(True, linestyle='--', linewidth=1, color='#555555')

        plt.show()
