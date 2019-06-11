# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
from urllib.parse import urlencode
from io import BytesIO

from astropy.utils.data import download_file
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import ICRS, Galactic, BaseCoordinateFrame
from astropy.coordinates import SkyCoord
from astropy import wcs

import cdshealpix
try:
    from astropy_healpix import HEALPix
except ImportError:
    pass

from ..abstract_moc import AbstractMOC
from ..interval_set import IntervalSet

from .. import core

from .boundaries import Boundaries
from .plot import fill, border

__author__ = "Thomas Boch, Matthieu Baumann"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, matthieu.baumann@astro.unistra.fr"

class MOC(AbstractMOC):
    """
    Multi-order spatial coverage class
    
    A MOC describes the coverage of an arbitrary region on the unit sphere.
    MOCs are usually used for describing the global coverage of catalog/image surveys such as GALEX or SDSS.
    A MOC corresponds to a list of `HEALPix <https://healpix.sourceforge.io/>`__ cells at different depths.
    This class gives you the possibility to:
    
    1. Define `~mocpy.moc.MOC` objects:

    - From a FITS file that stores HEALPix cells (see `from_fits`).
    - Directly from a list of HEALPix cells expressed either as a numpy structural array (see `from_healpix_cells`) or a simple
      python dictionnary (see `from_json`).
    - From a list of sky coordinates (see `from_skycoords`, `from_lonlat`).
    - From a convex/concave polygon (see `from_polygon`).
    - From a cone (will be implemented in a next version).

    2. Perform fast logical operations between `~mocpy.moc.MOC` objects:

    - The `intersection`
    - The `union`
    - The `difference`
    - The `complement`

    3. Plot the `~mocpy.moc.MOC` objects:

    - Draw the MOC with its HEALPix cells (see `fill`)
    - Draw the perimeter of a MOC (see `border`)

    4. Get the sky coordinates defining the border(s) of `~mocpy.moc.MOC` objects (see `get_boundaries`).

    5. Serialize `~mocpy.moc.MOC` objects to `astropy.io.fits.HDUList` or JSON dictionary and save it to a file.
    """

    def __init__(self, interval_set=None):
        super(MOC, self).__init__(interval_set)
        self._fits_header_keywords = {'COORDSYS': ('C', 'reference frame (C=ICRS)')}

    def contains(self, ra, dec, keep_inside=True):
        """
        Returns a boolean mask array of the positions lying inside (or outside) the MOC instance.

        Parameters
        ----------
        ra : `astropy.units.Quantity`
            Right ascension array
        dec : `astropy.units.Quantity`
            Declination array
        keep_inside : bool, optional
            True by default. If so the mask describes coordinates lying inside the MOC. If ``keep_inside``
            is false, contains will return the mask of the coordinates lying outside the MOC.

        Returns
        -------
        array : `~np.ndarray`
            A mask boolean array
        """
        max_depth = self.max_order
        m = np.zeros(3 << (2*(max_depth + 1)), dtype=bool)

        pix_id = core.flatten_pixels(self._interval_set._intervals, max_depth)
        m[pix_id] = True

        if not keep_inside:
            m = np.logical_not(m)

        pix = cdshealpix.lonlat_to_healpix(ra, dec, max_depth)

        return m[pix]

    def add_neighbours(self):
        """
        Extends the MOC instance so that it includes the HEALPix cells touching its border.

        The depth of the HEALPix cells added at the border is equal to the maximum depth of the MOC instance.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            self extended by one degree of neighbours.
        """
        max_depth = self.max_order

        # Get the pixels array of the MOC at the its max order.
        ipix = core.flatten_pixels(self._interval_set._intervals, max_depth)
        # Get the HEALPix array containing the neighbors of ``ipix``.
        # This array "extends" ``ipix`` by one degree of neighbors.
        ipix_extended = cdshealpix.neighbours(ipix, max_depth)
        ipix_extended = ipix_extended[ipix_extended >= 0]
        ipix_extended = ipix_extended.astype(np.uint64)

        # Compute the difference between ``extend_ipix`` and ``ipix`` to get only the neighboring pixels
        # located at the border of the MOC.
        ipix_neighbors = np.setdiff1d(ipix_extended, ipix)

        depth_neighbors = np.full(shape=ipix_neighbors.shape, fill_value=max_depth, dtype=np.int8)
        #intervals_neighbors = core.from_healpix_cells(ipix_neighbors, depth_neighbors)
        moc_neighbors = MOC.from_healpix_cells(ipix_neighbors, depth_neighbors)
        
        # This array of HEALPix neighbors are added to the MOC to get an ``extended`` MOC
        #self._interval_set._intervals = core.union(self._interval_set._intervals, moc_neighbors._interval_set._intervals)
        final = self.union(moc_neighbors)
        self._interval_set = final._interval_set
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
        max_depth = self.max_order
        # Get the HEALPix cells of the MOC at its max depth
        ipix = core.flatten_pixels(self._interval_set._intervals, max_depth)
        # Get the HEALPix array containing the neighbors of ``ipix``.
        # This array "extends" ``ipix`` by one degree of neighbors.
        ipix_extended = cdshealpix.neighbours(ipix, max_depth)
        ipix_extended = ipix_extended[ipix_extended >= 0]
        ipix_extended = ipix_extended.astype(np.uint64)
        # Get only the max depth HEALPix cells lying at the border of the MOC
        ipix_neighbors = np.setxor1d(ipix_extended, ipix)

        # Remove these pixels from ``ipix``
        ipix_around_border = cdshealpix.neighbours(ipix_neighbors, max_depth)
        ipix_around_border = ipix_around_border[ipix_around_border >= 0]
        ipix_around_border = ipix_around_border.astype(np.uint64)

        final_ipix = np.setdiff1d(ipix, ipix_around_border)
        final_depth = np.full(shape=final_ipix.shape, fill_value=max_depth, dtype=np.int8)

        # Build the reduced MOC, i.e. MOC without its pixels which were located at its border.
        final = MOC.from_healpix_cells(final_ipix, final_depth)
        self._interval_set = final._interval_set
        return self

    def fill(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Draws the MOC on a matplotlib axis.

        This performs the projection of the cells from the world coordinate system to the pixel image coordinate system.
        You are able to specify various styling kwargs for `matplotlib.patches.PathPatch`
        (see the `list of valid keywords <https://matplotlib.org/api/_as_gen/matplotlib.patches.PathPatch.html#matplotlib.patches.PathPatch>`__).

        Parameters
        ----------
        ax : `matplotlib.axes.Axes`
            Matplotlib axis.
        wcs : `astropy.wcs.WCS`
            WCS defining the World system <-> Image system projection.
        kw_mpl_pathpatch
            Plotting arguments for `matplotlib.patches.PathPatch`.
        
        Examples
        --------
        >>> from mocpy import MOC, WCS
        >>> from astropy.coordinates import Angle, SkyCoord
        >>> import astropy.units as u
        >>> # Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
        >>> filename = './../resources/P-GALEXGR6-AIS-FUV.fits'
        >>> moc = MOC.from_fits(filename)
        >>> # Plot the MOC using matplotlib
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure(111, figsize=(15, 15))
        >>> # Define a WCS as a context
        >>> with WCS(fig, 
        ...         fov=50 * u.deg,
        ...         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
        ...         coordsys="icrs",
        ...         rotation=Angle(0, u.degree),
        ...         projection="AIT") as wcs:
        ...     ax = fig.add_subplot(1, 1, 1, projection=wcs)
        ...     # Call fill giving the matplotlib axe and the `~astropy.wcs.WCS` object.
        ...     # We will set the matplotlib keyword linewidth to 0 so that it does not plot
        ...     # the border of each HEALPix cell.
        ...     # The color can also be specified along with an alpha value.
        ...     moc.fill(ax=ax, wcs=wcs, linewidth=0, alpha=0.5, fill=True, color="green")
        >>> plt.xlabel('ra')
        >>> plt.ylabel('dec')
        >>> plt.grid(color="black", linestyle="dotted")
        """
        fill.fill(self, ax, wcs, **kw_mpl_pathpatch)

    def border(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Draws the MOC border(s) on a matplotlib axis.

        This performs the projection of the sky coordinates defining the perimeter of the MOC to the pixel image coordinate system.
        You are able to specify various styling kwargs for `matplotlib.patches.PathPatch` 
        (see the `list of valid keywords <https://matplotlib.org/api/_as_gen/matplotlib.patches.PathPatch.html#matplotlib.patches.PathPatch>`__).

        Parameters
        ----------
        ax : `matplotlib.axes.Axes`
            Matplotlib axis.
        wcs : `astropy.wcs.WCS`
            WCS defining the World system <-> Image system projection. 
        kw_mpl_pathpatch
            Plotting arguments for `matplotlib.patches.PathPatch`

        Examples
        --------
        >>> from mocpy import MOC, WCS
        >>> from astropy.coordinates import Angle, SkyCoord
        >>> import astropy.units as u
        >>> # Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
        >>> filename = './../resources/P-GALEXGR6-AIS-FUV.fits'
        >>> moc = MOC.from_fits(filename)
        >>> # Plot the MOC using matplotlib
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure(111, figsize=(15, 15))
        >>> # Define a WCS as a context
        >>> with WCS(fig, 
        ...         fov=50 * u.deg,
        ...         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
        ...         coordsys="icrs",
        ...         rotation=Angle(0, u.degree),
        ...         projection="AIT") as wcs:
        ...     ax = fig.add_subplot(1, 1, 1, projection=wcs)
        ...     # Call border giving the matplotlib axe and the `~astropy.wcs.WCS` object.
        ...     moc.border(ax=ax, wcs=wcs, alpha=0.5, color="red")
        >>> plt.xlabel('ra')
        >>> plt.ylabel('dec')
        >>> plt.grid(color="black", linestyle="dotted")
        """
        border.border(self, ax, wcs, **kw_mpl_pathpatch)

    def get_boundaries(self, order=None):
        """
        Returns the sky coordinates defining the border(s) of the MOC.

        The border(s) are expressed as a list of SkyCoord.
        Each SkyCoord refers to the coordinates of one border of the MOC (i.e. 
        either a border of a connexe MOC part or a border of a hole
        located in a connexe MOC part).
        This function is currently not stable: encoding a vertice of a 
        HEALPix cell (N, E, S, W) should not depend on the position of the
        vertice but rather on the uniq value (+ 2 bits to encode the direction
        of the vertice).

        Parameters
        ----------
        order : int
            The depth of the MOC before computing its boundaries.
            A shallow depth leads to a faster computation.
            By default the maximum depth of the MOC is taken.

        Raises
        ------
        DeprecationWarning
            This method is not stable and not tested! A future more stable algorithm will be implemented!

        Return
        ------
        coords: [`~astropy.coordinates.SkyCoord`]
            A list of `~astropy.coordinates.SkyCoord` each describing one border.
        """

        import warnings
        warnings.simplefilter('default')
        warnings.warn('This method is not stable. A future more stable algorithm will be implemented!', DeprecationWarning)
        return Boundaries.get(self, order)

    @classmethod
    def from_image(cls, header, max_norder, mask=None):
        """
        Creates a `~mocpy.moc.MOC` from an image stored as a FITS file.

        Parameters
        ----------
        header : `astropy.io.fits.Header`
            FITS header containing all the info of where the image is located (position, size, etc...)
        max_norder : int
            The moc resolution.
        mask : `numpy.ndarray`, optional
            A boolean array of the same size of the image where pixels having the value 1 are part of
            the final MOC and pixels having the value 0 are not.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The resulting MOC.
        """
        # load the image data
        height = header['NAXIS2']
        width = header['NAXIS1']

        # use wcs from astropy to locate the image in the world coordinates
        w = wcs.WCS(header)

        if mask is not None:
            # We have an array of pixels that are part of of survey
            y, x = np.where(mask)
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
        world_crd = SkyCoord(w.wcs_pix2world(pix_crd, 1), unit="deg", frame=frame).icrs

        lon = world_crd.ra
        lat = world_crd.dec
        moc = MOC.from_lonlat(lon=lon, lat=lat, max_norder=max_norder)
        return moc

    @classmethod
    def from_fits_images(cls, path_l, max_norder):
        """
        Loads a MOC from a set of FITS file images.

        Parameters
        ----------
        path_l : [str]
            A list of path where the fits image are located.
        max_norder : int
            The MOC resolution.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The union of all the MOCs created from the paths found in ``path_l``.
        """
        moc = None
        for path in path_l:
            header = fits.getheader(path)
            current_moc = MOC.from_image(header=header, max_norder=max_norder)
            if moc is None:
                moc = current_moc
            else:
                moc = moc.union(current_moc)

        return moc

    @classmethod
    def from_vizier_table(cls, table_id, nside=256):
        """
        Creates a `~mocpy.moc.MOC` object from a VizieR table.

        **Info**: This method is already implemented in `astroquery.cds <https://astroquery.readthedocs.io/en/latest/cds/cds.html>`__. You can ask to get a `mocpy.moc.MOC` object
        from a vizier catalog ID.

        Parameters
        ----------
        table_id : str
            table index
        nside : int, optional
            256 by default

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC.
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
        Creates a `~mocpy.moc.MOC` object from a given ivorn.

        Parameters
        ----------
        ivorn : str
        nside : int, optional
            256 by default

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC.
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
        Creates a `~mocpy.moc.MOC` object from a given url.

        Parameters
        ----------
        url : str
            The url of a FITS file storing a MOC.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC.
        """
        path = download_file(url, show_progress=False, timeout=60)
        return cls.from_fits(path)

    @classmethod
    def from_skycoords(cls, skycoords, max_norder):
        """
        Creates a MOC from an `astropy.coordinates.SkyCoord`.

        Parameters
        ----------
        skycoords : `astropy.coordinates.SkyCoord`
            The sky coordinates that will belong to the MOC.
        max_norder : int
            The depth of the smallest HEALPix cells contained in the MOC.
        
        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        return cls.from_lonlat(lon=skycoords.icrs.ra, lat=skycoords.icrs.dec, max_norder=max_norder)

    @classmethod
    def from_lonlat(cls, lon, lat, max_norder):
        """
        Creates a MOC from astropy lon, lat `astropy.units.Quantity`.
        
        Parameters
        ----------
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates belonging to the MOC.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates belonging to the MOC.
        max_norder : int
            The depth of the smallest HEALPix cells contained in the MOC.
        
        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        intervals = core.from_lonlat(max_norder, lon.to_value(u.rad).astype(np.float64), lat.to_value(u.rad).astype(np.float64))
        return cls(IntervalSet(intervals, make_consistent=False))

    @classmethod
    def from_polygon_skycoord(cls, skycoord, max_depth=10):
        """
        Creates a MOC from a polygon.

        The polygon is given as an `astropy.coordinates.SkyCoord` that contains the 
        vertices of the polygon. Concave, convex and self-intersecting polygons are accepted.

        Parameters
        ----------
        skycoord : `astropy.coordinates.SkyCoord`
            The sky coordinates defining the vertices of a polygon. It can describe a convex or
            concave polygon but not a self-intersecting one.
        max_depth : int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        return MOC.from_polygon(lon=skycoord.icrs.ra, lat=skycoord.icrs.dec, max_depth=max_depth)

    @classmethod
    def from_polygon(cls, lon, lat, max_depth=10):
        """
        Creates a MOC from a polygon

        The polygon is given as lon and lat `astropy.units.Quantity` that define the 
        vertices of the polygon. Concave, convex and self-intersecting polygons are accepted.

        Parameters
        ----------
        lon : `astropy.units.Quantity`
            The longitudes defining the polygon. Can describe convex and
            concave polygons but not self-intersecting ones.
        lat : `astropy.units.Quantity`
            The latitudes defining the polygon. Can describe convex and concave
            polygons but not self-intersecting ones.
        max_depth : int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        pix, depth, fully_covered_flags = cdshealpix.polygon_search(lon, lat, max_depth)
        return MOC.from_healpix_cells(pix, depth, fully_covered_flags)

    @classmethod
    def from_healpix_cells(cls, ipix, depth, fully_covered=None):
        """
        Creates a MOC from a set of HEALPix cells at a given depth.

        Parameters
        ----------
        ipix : `numpy.ndarray`
            HEALPix cell indices. dtype must be np.uint64
        depth : `numpy.ndarray`
            Depth of the HEALPix cells. Must be of the same size of `ipix`.
            dtype must be np.uint8
        fully_covered : `numpy.ndarray`, optional
            HEALPix cells coverage flags. This flag informs whether a cell is
            fully covered by a cone (resp. polygon, elliptical cone) or not.
            Must be of the same size of `ipix`.

        Raises
        ------
        IndexError
            When `ipix`, `depth` and `fully_covered` do not have the same shape

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The MOC
        """
        if ipix.shape != depth.shape:
            raise IndexError("pixels and depth arrays must have the same shape")

        if fully_covered is not None and fully_covered.shape != ipix.shape:
            raise IndexError("fully covered and depth arrays must have the same shape")

        intervals = core.from_healpix_cells(ipix.astype(np.uint64), depth.astype(np.int8))
        return cls(IntervalSet(intervals, make_consistent=False))

    @property
    def sky_fraction(self):
        """
        Sky fraction covered by the MOC
        """
        max_depth = self.max_order
        flattened_pixels = core.flatten_pixels(self._interval_set._intervals, max_depth)
        num_pixels = flattened_pixels.size
        return num_pixels / float(3 << (2*(max_depth + 1)))

    # TODO : move this in astroquery.Simbad.query_region
    # See https://github.com/astropy/astroquery/pull/1466
    def query_simbad(self, max_rows=10000):
        """
        Query a view of SIMBAD data for SIMBAD objects in the coverage of the MOC instance.
        """
        return self._query('SIMBAD', max_rows)

    # TODO : move this in astroquery.Vizier.query_region
    # See https://github.com/astropy/astroquery/pull/1466
    def query_vizier_table(self, table_id, max_rows=10000):
        """
        Query a VizieR table for sources in the coverage of the MOC instance.
        """
        return self._query(table_id, max_rows)

    # TODO : move this in astroquery
    def _query(moc, resource_id, max_rows=100000):
        """
        Internal method to query Simbad or a VizieR table
        for sources in the coverage of the MOC instance
        """
        from astropy.io.votable import parse_single_table
        import requests

        moc_file = BytesIO()
        moc_fits = moc.serialize(format='fits')
        moc_fits.writeto(moc_file)

        r = requests.post('http://cdsxmatch.u-strasbg.fr/QueryCat/QueryCat',
                          data={'mode': 'mocfile',
                                'catName': resource_id,
                                'format': 'votable',
                                'limit': max_rows},
                          files={'moc': moc_file.getvalue()},
                          headers={'User-Agent': 'MOCPy'},
                          stream=True)

        votable = BytesIO()
        votable.write(r.content)

        table = parse_single_table(votable).to_table()

        return table

    def plot(self, title='MOC', frame=None):
        """
        Plot the MOC object using a mollweide projection.

        **Deprecated**: New `fill` and `border` methods produce more reliable results and allow you to specify additional 
        matplotlib style parameters.

        Parameters
        ----------
        title : str
            The title of the plot
        frame : `astropy.coordinates.BaseCoordinateFrame`, optional
            Describes the coordinate system the plot will be (ICRS, Galactic are the only coordinate systems supported).
        """
        import warnings
        warnings.simplefilter('default')
        warnings.warn('This method is deprecated and is no longer tested.'
                      'Please refer to this documentation page for plotting MOCs using'
                      'matplotlib: https://mocpy.readthedocs.io/en/latest/examples/examples.html#loading-and-plotting-the-moc-of-sdss', DeprecationWarning)

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

        m = np.zeros(12*4**(plotted_moc.max_order))
        pix_id = core.flatten_pixels(plotted_moc._interval_set._intervals, plotted_moc.max_order)

        # change the HEALPix cells if the frame of the MOC is not the same as the one associated with the plot method.
        if isinstance(frame, Galactic):
            lon, lat = hp.boundaries_lonlat(pix_id, step=2)
            sky_crd = SkyCoord(lon, lat, unit='deg')
            pix_id = hp.lonlat_to_healpix(sky_crd.galactic.l, sky_crd.galactic.b)

        m[pix_id] = 1

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
