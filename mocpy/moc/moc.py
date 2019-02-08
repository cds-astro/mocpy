# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from ..py23_compat import urlencode, int, BytesIO
import requests
import tempfile
import os
import numpy as np

from astropy.utils.data import download_file
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import ICRS, Galactic, BaseCoordinateFrame
from astropy.coordinates import SkyCoord
from astropy import wcs

from astropy_healpix import HEALPix
from astropy_healpix.healpy import nside2npix

from ..abstract_moc import AbstractMOC
from ..interval_set import IntervalSet


from .boundaries import Boundaries
from .plot import fill, border

__author__ = "Thomas Boch, Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, matthieu.baumann@astro.unistra.fr"

class MOC(AbstractMOC):
    """Multi-order spatial coverage class"""

    def __init__(self, interval_set=None):
        AbstractMOC.__init__(self, interval_set)
        self._fits_header_keywords = {'COORDSYS': ('C', 'reference frame (C=ICRS)')}

    def _best_res_pixels(self):
        """
        Get a numpy array of all the HEALPix indexes contained in the MOC at its max order.

        Returns
        -------
        array : `numpy.ndarray`
            the array of HEALPix at ``max_order``
        """
        factor = 2 * (AbstractMOC.HPY_MAX_NORDER - self.max_order)
        pix_l = []
        for iv in self._interval_set._intervals:
            for val in range(iv[0] >> factor, iv[1] >> factor):
                pix_l.append(val)

        return np.asarray(pix_l)

    def contains(self, ra, dec, keep_inside=True):
        """
        Get a boolean mask array of the positions lying inside (or outside) the MOC.

        Parameters
        ----------
        ra : `astropy.units.Quantity`
            right ascension array
        dec: `astropy.units.Quantity`
            declination array
        keep_inside : bool, optional
            True by default. If so the mask describes coordinates lying inside the MOC. If ``keep_inside``
            is false, the mask indicates the coordinates lying outside the MOC. 

        Returns
        -------
        array : `~numpy.darray`
            A mask boolean array
        """
        max_order = self.max_order
        m = np.zeros(nside2npix(1 << max_order), dtype=bool)

        pix_id_arr = self._best_res_pixels()
        m[pix_id_arr] = True

        if not keep_inside:
            m = np.logical_not(m)

        hp = HEALPix(nside=(1 << self.max_order), order='nested')
        pix_arr = hp.lonlat_to_healpix(ra, dec)

        return m[pix_arr]

    def _get_max_pix(self):
        return 3*(2**60)

    def add_neighbours(self):
        """
        Extends the MOC so that it includes the HEALPix cells touching the border of the MOC.

        This operation is done at the max depth of the MOC.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            self extended by one degree of neighbors
        """
        # Get the pixels array of the MOC at the its max order.
        ipix = self._best_res_pixels()

        hp = HEALPix(nside=(1 << self.max_order), order='nested')
        # Get the HEALPix array containing the neighbors of ``ipix``.
        # This array "extends" ``ipix`` by one degree of neighbors. 
        extend_ipix = AbstractMOC._neighbour_pixels(hp, ipix)
        
        # Compute the difference between ``extend_ipix`` and ``ipix`` to get only the neighboring pixels
        # located at the border of the MOC.
        neigh_ipix = np.setdiff1d(extend_ipix, ipix)

        shift = 2 * (AbstractMOC.HPY_MAX_NORDER - self.max_order)
        neigh_itv = np.vstack((neigh_ipix << shift, (neigh_ipix + 1) << shift)).T
        # This array of HEALPix neighbors are added to the MOC to get an ``extended`` MOC at its max order.
        self._interval_set = self._interval_set.union(IntervalSet.from_numpy_array(neigh_itv))
        return self

    def remove_neighbours(self):
        """
        Remove from the MOC the HEALPix cells located at the border of the MOC.

        This operation is done at the max depth of the MOC.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            self
        """
        # Get the HEALPix cells of the MOC at its max depth
        ipix = self._best_res_pixels()

        hp = HEALPix(nside=(1 << self.max_order), order='nested')
        # Extend it to include the max depth neighbor cells.
        extend_ipix = AbstractMOC._neighbour_pixels(hp, ipix)

        # Get only the max depth HEALPix cells lying at the border of the MOC
        neigh_ipix = np.setxor1d(extend_ipix, ipix)

        # Remove these pixels from ``ipix``
        border_ipix = AbstractMOC._neighbour_pixels(hp, neigh_ipix)
        reduced_ipix = np.setdiff1d(ipix, border_ipix)

        # Build the reduced MOC, i.e. MOC without its pixels which were located at its border.
        shift = 2 * (AbstractMOC.HPY_MAX_NORDER - self.max_order)
        reduced_itv = np.vstack((reduced_ipix << shift, (reduced_ipix + 1) << shift)).T
        self._interval_set = IntervalSet.from_numpy_array(reduced_itv)
        return self

    def fill(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Add the MOC to a matplotlib axis

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            Matplotlib axis
        wcs: `~astropy.wcs.WCS`
            World coordinate system in which the MOC is projeted
        kw_mpl_pathpatch
            Plotting arguments for `matplotlib.patches.PathPatch`
        """
        fill.fill(self, ax, wcs, **kw_mpl_pathpatch)

    def border(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Add the MOC border to a matplotlib axis

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            Matplotlib axis
        wcs: `~astropy.wcs.WCS`
            World coordinate system in which the MOC is projeted
        kw_mpl_pathpatch
            Plotting arguments for `matplotlib.patches.PathPatch`
        """
        border.border(self, ax, wcs, **kw_mpl_pathpatch)

    def get_boundaries(self, order=None):
        """
        Return the boundaries of the MOC

        The returned boundaries are expressed as a list of SkyCoord.
        Each SkyCoord refers to the coordinates of one border of the MOC (i.e. 
        either a border of a connexe MOC component or a border of a hole
        located in a connexe MOC component).

        Parameters
        ----------
        order: int
            The order of the MOC before computing its boundaries.
            A shallow order leads to a faster computation.
            By default the current max order of the MOC is taken.

        Return
        ------
        boundaries: [`~astropy.coordinates.SkyCoord`]
            A list of SkyCoords each describing one border.
        """
        return Boundaries.get(self, order)

    @classmethod
    def from_image(cls, header, max_norder, mask_arr=None):
        """
        Create a `~mocpy.moc.MOC` from an image stored as a fits file

        Parameters
        ----------
        header : `~astropy.io.fits.Header`
            fits header containing all the info of where the image is located (position, size, etc...)
        max_norder : int
            the moc resolution
        mask_arr : `~numpy.ndarray`, optional
            a 2D boolean array of the same size of the image where pixels having the value 1 are part of
            the final MOC and pixels having the value 0 are not.

        Returns
        -------
        moc : `mocpy.MOC`
            the MOC object loaded from the ``mask_arr`` and ``header`` extracted from the image
        """
        # load the image data
        height = header['NAXIS2']
        width = header['NAXIS1']

        # use wcs from astropy to locate the image in the world coordinates
        w = wcs.WCS(header)

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

        frame = wcs.utils.wcs_to_celestial_frame(w)
        world_pix_crd = SkyCoord(w.wcs_pix2world(pix_crd, 1), unit='deg', frame=frame)

        hp = HEALPix(nside=(1 << max_norder), order='nested', frame=ICRS())
        ipix = hp.skycoord_to_healpix(world_pix_crd)
        # remove doubles
        ipix = np.unique(ipix)

        shift = 2 * (AbstractMOC.HPY_MAX_NORDER - max_norder)
        intervals_arr = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        # This MOC will be consistent when one will do operations on the moc (union, inter, ...) or
        # simply write it to a fits or json file
        interval_set = IntervalSet.from_numpy_array(intervals_arr)

        return cls(interval_set=interval_set)

    @classmethod
    def from_fits_images(cls, path_l, max_norder):
        """
        Load a moc from a set of fits images

        Parameters
        ----------
        path_l : [str]
            the path list where the fits image are located
        max_norder : int
            moc resolution

        Returns
        -------
        moc : `~mocpy.MOC`
            the union of all the moc from path_l
        """
        moc = MOC()
        for path in path_l:
            header = fits.getheader(path)
            current_moc = MOC.from_image(header=header, max_norder=max_norder)
            moc = moc.union(current_moc)

        return moc

    @classmethod
    def from_vizier_table(cls, table_id, nside=256):
        """
        Create a `~mocpy.moc.MOC` object from a VizieR table

        TODO already implemented in astroquery.cds (asking the MOCServer for an ID).
        TODO a astroquery.cds.query_object(ID) should return a Dataset object containing the field ``moc_access_url``

        Parameters
        ----------
        table_id : str
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

        TODO already implemented in astroquery.cds (asking the MOCServer for an ID).
        TODO a astroquery.cds.query_object(ID) should return a Dataset object containing the field ``moc_access_url``

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
        path = download_file(url, show_progress=False, timeout=60)
        return cls.from_fits(path)

    @classmethod
    def from_skycoords(cls, skycoords, max_norder):
        """
        Create a MOC from an astropy skycoord.

        Parameters
        ----------
        skycoords : `~astropy.coordinates.SkyCoord`
            Must be expressed using the ICRS coordinate system.
        max_norder : int
            The depth of the smallest HEALPix cells contained in the MOC.
        
        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        hp = HEALPix(nside=(1 << max_norder), order='nested')
        ipix = hp.lonlat_to_healpix(skycoords.icrs.ra, skycoords.icrs.dec)

        shift = 2 * (AbstractMOC.HPY_MAX_NORDER - max_norder)
        intervals_arr = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        interval_set = IntervalSet.from_numpy_array(intervals_arr)
        return cls(interval_set)

    @classmethod
    def from_lonlat(cls, lon, lat, max_norder):
        """
        Create a MOC from astropy lon, lat quantities.
        
        Parameters
        ----------
        skycoords : `~astropy.coordinates.SkyCoord`
            Must be expressed in ICRS because the astropy-HEALPix frame used here
            is in ICRS.
        max_norder : int
            The depth of the smallest HEALPix cells contained in the MOC.
        
        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        hp = HEALPix(nside=(1 << max_norder), order='nested')
        ipix = hp.lonlat_to_healpix(lon, lat)

        shift = 2 * (AbstractMOC.HPY_MAX_NORDER - max_norder)
        intervals_arr = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        interval_set = IntervalSet.from_numpy_array(intervals_arr)
        return cls(interval_set)

    @classmethod
    def from_polygon_skycoord(cls, skycoord, inside=None, max_depth=10):
        """
        Create a MOC from a polygon

        Parameters
        ----------
        skycoord : `~astropy.coordinates.SkyCoord`
            The sky coordinates defining the polygon. Can describe convex and
            concave polygons but not self-intersecting ones.
        inside : `~astropy.coordinates.SkyCoord`, optional
            A point that will be inside the MOC is needed as it is not possible to determine the inside area of a polygon 
            on the unit sphere (there is no infinite area that can be considered as the outside.
            On the sphere, a polygon delimits two finite areas.).
            Possible improvement: take the inside area as the one covering the less of the sphere.

            If inside=None (default behavior), the mean of all the vertices is taken as lying inside the polygon. That approach may not work for 
            concave polygons.
        max_depth: int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        return MOC.from_polygon(lon=skycoord.icrs.ra, lat=skycoord.icrs.dec,
                                inside=inside, max_depth=max_depth)

    @classmethod
    def from_polygon(cls, lon, lat, inside=None, max_depth=10):
        """
        Create a MOC from a polygon

        Parameters
        ----------
        lon : `~astropy.units.Quantity`
            The longitudes defining the polygon. Can describe convex and
            concave polygons but not self-intersecting ones.
        lat : `~astropy.units.Quantity`
            The latitudes defining the polygon. Can describe convex and concave
            polygons but not self-intersecting ones.
        inside : `~astropy.coordinates.SkyCoord`, optional
            A point that will be inside the MOC is needed as it is not possible to determine the inside area of a polygon 
            on the unit sphere (there is no infinite area that can be considered as the outside.
            On the sphere, a polygon delimits two finite areas.).
            Possible improvement: take the inside area as the one covering the less of the sphere.

            If inside=None (default behavior), the mean of all the vertices is taken as lying inside the polygon. That approach may not work for 
            concave polygons.
        max_depth: int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        from .polygon import PolygonComputer

        polygon_computer = PolygonComputer(lon, lat, inside, max_depth)
        # Create the moc from the python dictionary

        moc = MOC.from_json(polygon_computer.ipix)
        # We degrade it to the user-requested order
        if polygon_computer.degrade_to_max_depth:
            moc = moc.degrade_to_order(max_depth)

        return moc

    @property
    def sky_fraction(self):
        """
        Return the sky fraction covered by the MOC
        """
        pix_id_arr = self._best_res_pixels()
        nb_pix_filled = pix_id_arr.size
        return nb_pix_filled / float(3 << (2*(self.max_order + 1)))

    # TODO : move this in astroquery.Simbad.query_region
    def query_simbad(self, max_rows=10000):
        """
        Query a view of SIMBAD data for SIMBAD objects in the coverage of the MOC instance
        """
        return self._query('SIMBAD', max_rows)

    # TODO : move this in astroquery.Vizier.query_region
    def query_vizier_table(self, table_id, max_rows=10000):
        """
        Query a VizieR table for sources in the coverage of the MOC instance
        """
        return self._query(table_id, max_rows)

    # TODO : move this in astroquery
    def _query(self, resource_id, max_rows):
        """
        Internal method to query Simbad or a VizieR table
        for sources in the coverage of the MOC instance
        """
        from astropy.io.votable import parse_single_table

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

        tmp_vot = BytesIO()
        tmp_vot.write(r.content)

        table = parse_single_table(tmp_vot).to_table()

        # finally delete temp files
        os.unlink(tmp_moc.name)

        return table

    def plot(self, title='MOC', frame=None):
        """
        Plot the MOC object using a mollweide projection.

        This method uses matplotlib. New `fill` and `border` methods tend to replace
        this one, produce more reliable results and allow you to specify additional 
        matplotlib style rendering parameters.

        Parameters
        ----------
        title : str
            the title of the plot
        frame : `astropy.coordinates.BaseCoordinateFrame`, optional
            Describes the coordinate system the plot will be (ICRS, Galactic are the only coordinate systems supported).
        """
        frame = ICRS() if frame is None else frame

        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.pyplot as plt

        plot_order = 8
        if self.max_order > plot_order:
            plotted_moc = self.degrade_to_order(plot_order)
        else:
            plotted_moc = self

        num_pixels_map = 1024
        delta = 2. * np.pi / num_pixels_map

        x = np.arange(-np.pi, np.pi, delta)
        y = np.arange(-np.pi/2, np.pi/2, delta)
        lon_rad, lat_rad = np.meshgrid(x, y)
        hp = HEALPix(nside=(1 << plotted_moc.max_order), order='nested')

        if frame and not isinstance(frame, BaseCoordinateFrame):
            raise ValueError("Only Galactic/ICRS coordinate systems are supported."
                             "Please set `coord` to either 'C' or 'G'.")

        pix_map = hp.lonlat_to_healpix(lon_rad * u.rad, lat_rad * u.rad)

        m = np.zeros(nside2npix(1 << plotted_moc.max_order))
        pix_id_arr = plotted_moc._best_res_pixels()

        # change the HEALPix cells if the frame of the MOC is not the same as the one associated with the plot method.
        if isinstance(frame, Galactic):
            lon, lat = hp.boundaries_lonlat(pix_id_arr, step=2)
            sky_crd = SkyCoord(lon, lat, unit='deg')
            pix_id_arr = hp.lonlat_to_healpix(sky_crd.galactic.l, sky_crd.galactic.b)

        m[pix_id_arr] = 1

        z = np.flip(m[pix_map], axis=1)

        plt.figure(figsize=(10, 10))

        ax = plt.subplot(111, projection="mollweide")
        ax.set_xticklabels(['150°', '120°', '90°', '60°', '30°', '0°', '330°', '300°', '270°', '240°', '210°', '180°'])

        color_map = LinearSegmentedColormap.from_list('w2r', ['#eeeeee', '#aa0000'])
        color_map.set_under('w')
        color_map.set_bad('gray')

        ax.pcolormesh(x, y, z, cmap=color_map, vmin=0, vmax=1)
        ax.tick_params(labelsize=14, labelcolor='#000000')
        plt.title(title)
        plt.grid(True, linestyle='--', linewidth=1, color='#555555')

        plt.show()
