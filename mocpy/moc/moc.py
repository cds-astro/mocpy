# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
from urllib.parse import urlencode
from io import BytesIO

from astropy.utils.data import download_file
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import ICRS, Galactic, BaseCoordinateFrame
from astropy.coordinates import SkyCoord, Angle, Longitude, Latitude
from astropy import wcs

import cdshealpix
try:
    from astropy_healpix import HEALPix
except ImportError:
    pass

from ..abstract_moc import AbstractMOC
from ..interval_set import IntervalSet

from .. import mocpy

from .boundaries import Boundaries
from .plot import fill, border

__author__ = "Thomas Boch, Matthieu Baumann, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "thomas.boch@astro.unistra.fr, baumannmatthieu0@gmail.com, francois-xavier.pineau@astro.unistra.fr"

class MOC(AbstractMOC):
    """
    Multi-order spatial coverage class.

    A MOC describes the coverage of an arbitrary region on the unit sphere.
    MOCs are usually used for describing the global coverage of catalog/image surveys such as GALEX or SDSS.
    A MOC corresponds to a list of `HEALPix <https://healpix.sourceforge.io/>`__ cells at different depths.
    This class gives you the possibility to:
    
    1. Define `~mocpy.moc.MOC` objects:

    - From a FITS file that stores HEALPix cells (see `load(path, 'fits')`).
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

    # I introduced, but do not like, the double `make_consistent` (MOC + IntervalSet)
    # but `coverage_merge_time_intervals` is no more genric
    # and I can't remove `make_consistent` from `IntervalSet` without changing tests
    def __init__(self, interval_set=None, make_consistent=True, min_depth=None):
        """
        Moc constructor.

        The merging step of the overlapping intervals is done here.

        Parameters
        ----------
        intervals : `~numpy.ndarray`
            a N x 2 numpy array representing the set of intervals.
        make_consistent : bool, optional
            True by default. Remove the overlapping intervals that makes
            a valid MOC (i.e. can be plot, serialized, manipulated).
        """
        super(MOC, self).__init__(interval_set)

        if make_consistent:
            if min_depth is None:
                min_depth = -1

            min_depth = np.int8(min_depth)
            self._merge_intervals(min_depth)

    def _merge_intervals(self, min_depth):
        if not self.empty():
            self._interval_set._intervals = mocpy.coverage_merge_hpx_intervals(self._interval_set._intervals, min_depth)

    @property
    def max_order(self):
        """
        Depth of the smallest HEALPix cells found in the MOC instance.
        """
        depth = mocpy.hpx_coverage_depth(self._interval_set._intervals)
        depth = np.uint8(depth)
        return depth

    def refine_to_order(self, min_depth):
        intervals = mocpy.coverage_merge_hpx_intervals(self._interval_set._intervals, min_depth)
        interval_set = IntervalSet(intervals, make_consistent=False)
        return MOC(interval_set, make_consistent=False)

    def complement(self):
        """
        Returns the complement of the MOC instance.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC.
        """
        intervals = mocpy.hpx_coverage_complement(self._interval_set._intervals)
        interval_set = IntervalSet(intervals, make_consistent=False)
        return MOC(interval_set, make_consistent=False)

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
        intervals = mocpy.hpx_coverage_expand(self.max_order, self._interval_set._intervals)
        interval_set = IntervalSet(intervals, make_consistent=False)
        return MOC(interval_set, make_consistent=False)

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
        intervals = mocpy.hpx_coverage_contract(self.max_order, self._interval_set._intervals)
        interval_set = IntervalSet(intervals, make_consistent=False)
        return MOC(interval_set, make_consistent=False)

    def split_count(self, include_indirect_neighbours=False):
       """
       Returns the number of disjoint MOCs the given MOC contains.

       Parameter
       ---------
       include_indirect_neighbours : bool
           if `false`, only consider  cells having a common edge as been part of a same MOC
           if `true`, also consider cells having a common vertex as been part of the same MOC
       """
       return mocpy.hpx_coverage_split_count(self.max_order, include_indirect_neighbours, self._interval_set._intervals)

    def split(self, include_indirect_neighbours=False):
        """
        Returns the disjoint MOCs this MOC contains.format

        Parameter
        ---------
        include_indirect_neighbours : bool
            if `false`, only consider  cells having a common edge as been part of a same MOC
            if `true`, also consider cells having a common vertex as been part of the same MOC

        WARNING
        -------
        Please use `~mocpy.moc.MOC.split_count` first to ensure the number is not too high
        
        """
        list_of_intervals = mocpy.hpx_coverage_split(self.max_order, include_indirect_neighbours, self._interval_set._intervals)
        mocs = map(lambda intervals: MOC(IntervalSet(intervals, make_consistent=False), make_consistent=False), list_of_intervals)
        return mocs

    def degrade_to_order(self, new_order):
        """
        Degrades the MOC instance to a new, less precise, MOC.

        The maximum depth (i.e. the depth of the smallest HEALPix cells that can be found in the MOC) of the
        degraded MOC is set to ``new_order``.

        Parameters
        ----------
        new_order : int
            Maximum depth of the output degraded MOC.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The degraded MOC.
        """
        intervals = mocpy.hpx_coverage_degrade(self._interval_set._intervals, new_order)
        return MOC(IntervalSet(intervals, make_consistent=False), make_consistent=False)

    def contains(self, ra, dec, keep_inside=True):
        """
        Returns a boolean mask array of the positions lying inside (or outside) the MOC instance.

        Parameters
        ----------
        ra : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            Right ascension array
        dec : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
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

        pix_id = mocpy.flatten_pixels(self._interval_set._intervals, max_depth)
        m[pix_id] = True

        if not keep_inside:
            m = np.logical_not(m)

        ra  = ra  if isinstance(ra, Longitude) else Longitude(ra)
        dec = dec if isinstance(dec, Latitude) else Latitude(dec)

        pix = cdshealpix.lonlat_to_healpix(ra, dec, max_depth)

        return m[pix]

    ## TODO: implement: def contains_including_surrounding(self, ra, dec, distance)

    def add_neighbours(self):
        """
        Extends the MOC instance so that it includes the HEALPix cells touching its border.

        The depth of the HEALPix cells added at the border is equal to the maximum depth of the MOC instance.

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            self extended by one degree of neighbours.
        """
        intervals = mocpy.hpx_coverage_expand(self.max_order, self._interval_set._intervals)
        self._interval_set = IntervalSet(intervals, make_consistent=False)
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
        intervals = mocpy.hpx_coverage_contract(self.max_order, self._interval_set._intervals);
        self._interval_set = IntervalSet(intervals, make_consistent=False)
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
        >>> from mocpy import MOC, World2ScreenMPL
        >>> from astropy.coordinates import Angle, SkyCoord
        >>> import astropy.units as u
        >>> # Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
        >>> filename = './../resources/P-GALEXGR6-AIS-FUV.fits'
        >>> moc = MOC.load(filename, 'fits')
        >>> # Plot the MOC using matplotlib
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure(111, figsize=(15, 15))
        >>> # Define a WCS as a context
        >>> with World2ScreenMPL(fig, 
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
        >>> from mocpy import MOC, World2ScreenMPL
        >>> from astropy.coordinates import Angle, SkyCoord
        >>> import astropy.units as u
        >>> # Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
        >>> filename = './../resources/P-GALEXGR6-AIS-FUV.fits'
        >>> moc = MOC.load(filename, 'fits')
        >>> # Plot the MOC using matplotlib
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure(111, figsize=(15, 15))
        >>> # Define a WCS as a context
        >>> with World2ScreenMPL(fig, 
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

        Returns
        -------
        coords: [`~astropy.coordinates.SkyCoord`]
            A list of `~astropy.coordinates.SkyCoord` each describing one border.
        """
        import warnings
        warnings.warn('This method is not stable. A future more stable algorithm will be implemented!', DeprecationWarning)
        return Boundaries.get(self, order)

    @classmethod
    def from_fits_image(cls, hdu, max_norder, mask=None):
        """
        Creates a `~mocpy.moc.MOC` from an image stored as a FITS file.

        Parameters
        ----------
        hdu : HDU object
            HDU containing the data of the image
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
        # Only take the first HDU
        header = hdu.header

        height = header['NAXIS2']
        width = header['NAXIS1']

        # Compute a WCS from the header of the image
        w = wcs.WCS(header)

        if mask is None:
            data = hdu.data
            # A mask is computed discarding nan floating values
            mask = np.isfinite(data)

            # If the BLANK keyword is set to a value then we mask those
            # pixels too
            if header.get('BLANK') is not None:
                discard_val = header['BLANK']

                # We keep the finite values and those who are not equal to the BLANK field
                mask = mask & (data != discard_val)

        y, x = np.where(mask)
        pix = np.dstack((x, y))[0]

        world = w.wcs_pix2world(pix, 0)

        # Remove coord containing inf/nan values
        good = np.isfinite(world)

        # It is a good coordinates whether both its coordinate are good
        good = good[:, 0] & good[:, 1]
        world = world[good]

        # Get the frame from the wcs
        frame = wcs.utils.wcs_to_celestial_frame(w)
        skycrd = SkyCoord(
            world,
            unit="deg",
            frame=frame
        )

        # Compute the order based on the CDELT
        c1 = header['CDELT1']
        c2 = header['CDELT2']
        max_res_px = np.sqrt(c1*c1 + c2*c2) * np.pi / 180.0
        max_depth_px = int(np.floor(np.log2(np.pi / (3 * max_res_px * max_res_px)) / 2))

        max_norder = min(max_norder, max_depth_px)

        moc = MOC.from_lonlat(
            lon=skycrd.icrs.ra,
            lat=skycrd.icrs.dec,
            max_norder=max_norder
        )
        return moc

    @classmethod
    def from_fits_images(cls, path_l, max_norder):
        """
        Loads a MOC from a set of FITS file images.

        Assumes the data of the image is stored in the first HDU of the FITS file.
        Please call `~mocpy.moc.MOC.from_fits_image` for passing another hdu than the first one.

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
        moc = MOC()
        for filename in path_l:
            with fits.open(filename) as hdul:
                current_moc = MOC.from_fits_image(hdu=hdul[0], max_norder=max_norder)
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
        return cls.load(path, 'fits')

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
        lon : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            The longitudes of the sky coordinates belonging to the MOC.
        lat : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
            The latitudes of the sky coordinates belonging to the MOC.
        max_norder : int
            The depth of the smallest HEALPix cells contained in the MOC.
        
        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        intervals = mocpy.from_lonlat(max_norder, lon.to_value(u.rad).astype(np.float64), lat.to_value(u.rad).astype(np.float64))
        return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)

    @classmethod
    def from_multiordermap_fits_file(cls, path,
        cumul_from=0.0, cumul_to=1.0, 
        asc=False, strict=True, no_split=True, reverse_decent=False):
        """
        Creates a MOC from a mutli-order map FITS file.

        HEALPix cells are first sorted by their values.
        The MOC contains the cells from which the cumulative value is between
        ``cumul_from`` and ``cumul_to``.
        Cells being on the fence are recursively splitted and added
        until the depth of the cells is equal to ``max_norder``.

        For compatibility with Aladin, use ``no_split=False`` and ``reverse_decent=True``
        
        Remark: using ``no_split=False``, the way the cells overlapping with the low and high thresholds are split
        is somewhat arbitrary.

        Parameters
        ----------
        path : str
            The path to the file to save the MOC in.
        cumul_from : float
            Cumulative value from which cells will be added to the MOC
        cumul_to : float
            Cumulative value to which cells will be added to the MOC
        asc: boolean
            the cumulative value is computed from lower to highest densities instead of from highest to lowest
        strict: boolean
            (sub-)cells overlapping the `cumul_from` or `cumul_to` values are not added
        no_split: boolean
            cells overlapping the `cumul_from` or `cumul_to` values are not recursively split
        reverse_decent: boolean
            perform the recursive decent from the highest cell number to the lowest (to be compatible with Aladin)

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        intervals = mocpy.spatial_moc_from_multiordermap_fits_file(
            path,
            np.float64(cumul_from),
            np.float64(cumul_to),
            asc, 
            strict, 
            no_split, 
            reverse_decent
        )

        return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)

    @classmethod
    def from_valued_healpix_cells(cls, 
        uniq, values, max_depth=None, 
        cumul_from=0.0, cumul_to=1.0, 
        asc=False, strict=True, no_split=True, reverse_decent=False):
        """
        Creates a MOC from a list of uniq associated with values.

        HEALPix cells are first sorted by their values.
        The MOC contains the cells from which the cumulative value is between
        ``cumul_from`` and ``cumul_to``.
        Cells being on the fence are recursively splitted and added
        until the depth of the cells is equal to ``max_norder``.

        For compatibility with Aladin, use ``no_split=False`` and ``reverse_decent=True``
        
        Remark: using ``no_split=False``, the way the cells overlapping with the low and high thresholds are split
        is somewhat arbitrary.

        Parameters
        ----------
        uniq : `numpy.ndarray`
            HEALPix cell indices written in uniq. dtype must be np.uint64
        values : `numpy.ndarray`
            Probabilities associated with each ``uniq`` cells. dtype must be np.float64
        max_depth : int, optional
            The max depth of the MOC. If a depth is given, degrade the MOC to this depth before returning it to the user.
            Otherwise choose as ``max_depth`` the depth corresponding to the smallest HEALPix cell found in ``uniq``.
        cumul_from : float
            Cumulative value from which cells will be added to the MOC
        cumul_to : float
            Cumulative value to which cells will be added to the MOC
        asc: boolean
            the cumulative value is computed from lower to highest densities instead of from highest to lowest
        strict: boolean
            (sub-)cells overlapping the `cumul_from` or `cumul_to` values are not added
        no_split: boolean
            cells overlapping the `cumul_from` or `cumul_to` values are not recursively split
        reverse_decent: boolean
            perform the recursive decent from the highest cell number to the lowest (to be compatible with Aladin)

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        max_depth_tile = 0
        if uniq.size > 0:
            # Get the depth of the smallest uniq
            # Bigger uniq corresponds to big depth HEALPix cells.
            max_depth_tile = int(np.log2(uniq.max() >> 2)) >> 1
            assert max_depth_tile >= 0 and max_depth_tile <= 29, "Invalid uniq numbers. Too big uniq or negative uniq numbers might the cause."

        # Create the MOC at the max_depth equals to the smallest cell 
        # found in the uniq array
        intervals = mocpy.from_valued_hpx_cells(
            np.uint8(max_depth_tile),
            uniq.astype(np.uint64),
            values.astype(np.float64),
            np.float64(cumul_from),
            np.float64(cumul_to),
            asc, 
            strict, 
            no_split, 
            reverse_decent
        )
        moc = cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)

        # Degrade the MOC to the depth requested by the user
        if max_depth is not None:
            assert max_depth >= 0 and max_depth <= 29, "Max depth must be in [0, 29]"
            moc = moc.degrade_to_order(max_depth)

        return moc

    @classmethod
    def from_elliptical_cone(cls, lon, lat, a, b, pa, max_depth, delta_depth=2):
        """
        Creates a MOC from an elliptical cone

        The ellipse is centered around the (`lon`, `lat`) position. `a` (resp. `b`) corresponds
        to the semi-major axis magnitude (resp. semi-minor axis magnitude). `pa` is expressed as a
        `~astropy.coordinates.Angle` and defines the position angle of the elliptical cone.

        Parameters
        ----------
        lon : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            The longitude of the center of the elliptical cone.
        lat : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
            The latitude of the center of the elliptical cone.
        a : `astropy.coordinates.Angle`
            The semi-major axis angle of the elliptical cone.
        b : `astropy.coordinates.Angle`
            The semi-minor axis angle of the elliptical cone.
        pa : `astropy.coordinates.Angle`
            The position angle (i.e. the angle between the north and the semi-major axis, east-of-north).
        max_depth : int
            Maximum HEALPix cell resolution.
        delta_depth : int, optional
            To control the approximation, you can choose to perform the computations at a deeper
            depth using the `depth_delta` parameter.
            The depth at which the computations will be made will therefore be equal to
            `depth` + `depth_delta`.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC

        Examples
        --------
        >>> from mocpy import MOC
        >>> import astropy.units as u
        >>> from astropy.coordinates import Angle, Longitude, Latitude
        >>> moc = MOC.from_elliptical_cone(
        ...  lon=Longitude(0 * u.deg),
        ...  lat=Latitude(0 * u.deg),
        ...  a=Angle(10, u.deg),
        ...  b=Angle(5, u.deg),
        ...  pa=Angle(0, u.deg),
        ...  max_depth=10
        ... )
        """
        lon = lon if isinstance(lon, Longitude) else Longitude(lon) 
        lat = lat if isinstance(lat, Latitude)  else Latitude(lat) 
        pix, depth, fully_covered_flags = cdshealpix.elliptical_cone_search(lon, lat, a, b, pa, max_depth, delta_depth, flat=False)
        return MOC.from_healpix_cells(pix, depth, fully_covered_flags)

    @classmethod
    def from_cone(cls, lon, lat, radius, max_depth, delta_depth=2):
        """
        Creates a MOC from a cone.

        The cone is centered around the (`lon`, `lat`) position with a radius expressed by
        `radius`.

        Parameters
        ----------
        lon : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            The longitude of the center of the cone.
        lat : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
            The latitude of the center of the cone.
        radius : `astropy.coordinates.Angle`
            The radius angle of the cone.
        max_depth : int
            Maximum HEALPix cell resolution.
        delta_depth : int, optional
            To control the approximation, you can choose to perform the computations at a deeper
            depth using the `depth_delta` parameter.
            The depth at which the computations will be made will therefore be equal to
            `max_depth` + `depth_delta`.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC

        Examples
        --------
        >>> from mocpy import MOC
        >>> import astropy.units as u
        >>> from astropy.coordinates import Angle, Longitude, Latitude
        >>> moc = MOC.from_cone(
        ...  lon=Longitude(0 * u.deg),
        ...  lat=Latitude(0 * u.deg),
        ...  radius=Angle(10, u.deg),
        ...  max_depth=10
        ... )
        """
        lon = lon if isinstance(lon, Longitude) else Longitude(lon)
        lat = lat if isinstance(lat, Latitude)  else Latitude(lat)
        pix, depth, fully_covered_flags = cdshealpix.cone_search(lon, lat, radius, max_depth, delta_depth, flat=False)
        return MOC.from_healpix_cells(pix, depth, fully_covered_flags)

    @classmethod
    def from_ring(cls, lon, lat, internal_radius, external_radius, max_depth, delta_depth=2):
        """
        Creates a MOC from a ring.

        The cone is centered around the (`lon`, `lat`) position with an internal radius expressed by
        `internal_radius` and an external radius expressed by `external_radius`.

        Parameters
        ----------
        lon : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            The longitude of the center of the ring.
        lat : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
            The latitude of the center of the ring.
        internal_radius : `astropy.coordinates.Angle`
            The internal radius angle of the ring.
        external_radius : `astropy.coordinates.Angle`
            The external radius angle of the ring.
        max_depth : int
            Maximum HEALPix cell resolution.
        delta_depth : int, optional
            To control the approximation, you can choose to perform the computations at a deeper
            depth using the `depth_delta` parameter.
            The depth at which the computations will be made will therefore be equal to
            `max_depth` + `depth_delta`.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC

        Examples
        --------
        >>> from mocpy import MOC
        >>> import astropy.units as u
        >>> from astropy.coordinates import Angle, Longitude, Latitude
        >>> moc = MOC.from_ring(
        ...  lon=Longitude(0 * u.deg),
        ...  lat=Latitude(0 * u.deg),
        ...  internal_radius=Angle(5, u.deg),
        ...  external_radius=Angle(10, u.deg),
        ...  max_depth=10
        ... )
        """
        lon = lon if isinstance(lon, Longitude) else Longitude(lon)
        lat = lat if isinstance(lat, Latitude)  else Latitude(lat)
        intervals = mocpy.from_ring(
            np.float64(lon.degree),
            np.float64(lat.degree),
            np.float64(internal_radius.degree),
            np.float64(external_radius.degree),
            np.uint8(max_depth),
            np.uint8(delta_depth)
        )
        return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)

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
        lon : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            The longitudes defining the polygon. Can describe convex and
            concave polygons but not self-intersecting ones.
        lat : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
            The latitudes defining the polygon. Can describe convex and concave
            polygons but not self-intersecting ones.
        max_depth : int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        lon = lon if isinstance(lon, Longitude) else Longitude(lon)    
        lat = lat if isinstance(lat, Latitude)  else Latitude(lat)
        pix, depth, fully_covered_flags = cdshealpix.polygon_search(lon, lat, max_depth)
        return MOC.from_healpix_cells(pix, depth, fully_covered_flags)

    @classmethod
    def from_healpix_cells(cls, ipix, depth, fully_covered=None):
        """
        Creates a MOC from a set of HEALPix cells at a given depth.

        Parameters
        ----------
        ipix : `numpy.ndarray`
            HEALPix cell indices in the NESTED notation. dtype must be np.uint64
        depth : `numpy.ndarray`
            Depth of the HEALPix cells. Must be of the same size of `ipix`.
            dtype must be np.uint8. Corresponds to the `level` of an HEALPix cell in astropy.healpix.
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

        intervals = mocpy.from_healpix_cells(ipix.astype(np.uint64), depth.astype(np.uint8))
        return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)

    @staticmethod
    def order_to_spatial_resolution(order):
        """
        Convert a depth to its equivalent spatial resolution.

        Parameters
        ----------
        order : int
            Spatial depth.

        Returns
        -------
        spatial_resolution : `~astropy.coordinates.Angle`
            Spatial resolution.

        """
        spatial_resolution = Angle(np.sqrt(np.pi/(3 * 4**(order))), unit='rad')
        return spatial_resolution

    @staticmethod
    def spatial_resolution_to_order(spatial_resolution):
        """
        Convert a spatial resolution to a MOC order.

        Parameters
        ----------
        spatial_resolution : `~astropy.coordinates.Angle`
            Spatial resolution

        Returns
        -------
        order : int
            The order corresponding to the spatial resolution
        """
        res_rad = spatial_resolution.rad
        order = np.ceil(np.log2(np.pi/(3*res_rad*res_rad))/2)
        return np.uint8(order)

    @property
    def _fits_header_keywords(self):
        return {
            'PIXTYPE': 'HEALPIX',
            'ORDERING': 'NUNIQ',
            'COORDSYS': ('C', 'reference frame (C=ICRS)'),
            'MOCORDER': self.max_order,
            'MOCTOOL': 'MOCPy'
        }
    
    @property
    def _fits_format(self):
        depth = self.max_order
        if depth <= 13:
            fits_format = '1J'
        else:
            fits_format = '1K'
        return fits_format

    @property
    def sky_fraction(self):
        """
        Sky fraction covered by the MOC
        """
        sky_fraction = mocpy.coverage_sky_fraction(self._interval_set._intervals)
        return sky_fraction

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
        warnings.warn('This method is deprecated and is no longer tested.'
                      'Please refer to this documentation page for plotting MOCs using'
                      'matplotlib: https://cds-astro.github.io/mocpy/xamples/examples.html#loading-and-plotting-the-moc-of-sdss', DeprecationWarning)

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
        pix_id = mocpy.flatten_pixels(plotted_moc._interval_set._intervals, plotted_moc.max_order)

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

    def save(self, path, format='fits', overwrite=False):
        """
        Writes the Spatial MOC to a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str
            The path to the file to save the MOC in.
        format : str, optional
            The format in which the MOC is saved.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        overwrite : bool, optional
            If the file already exists and you want to overwrite it, then set the  ``overwrite`` keyword. 
            Default to False.
        """
        import os
        file_exists = os.path.isfile(path)
        
        if file_exists and not overwrite:
            raise OSError('File {} already exists! Set ``overwrite`` to '
                          'True if you want to replace it.'.format(path))
        
        if format == 'fits':
            mocpy.spatial_moc_to_fits_file(self.max_order, self._interval_set._intervals, path)
        elif format == 'ascii':
            mocpy.spatial_moc_to_ascii_file(self.max_order, self._interval_set._intervals, path)
        elif format == 'json':
            mocpy.spatial_moc_to_json_file(self.max_order, self._interval_set._intervals, path)
        else:
            formats = ('fits', 'ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))

    @classmethod
    def load(cls, path, format='fits'):
        """
        Load the Spatial MOC from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str
            The path to the file to load the MOC from.
        format : str, optional
            The format from which the MOC is loaded.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        """
        if format == 'fits':
            intervals = mocpy.spatial_moc_from_fits_file(path)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        elif format == 'ascii':
            intervals = mocpy.spatial_moc_from_ascii_file(path)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        elif format == 'json':
            intervals = mocpy.spatial_moc_from_json_file(path)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        else:
            formats = ('fits', 'ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))

    def to_string(self, format='ascii'):
        """
        Writes the Spatial MOC into a string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == 'ascii':
            return mocpy.spatial_moc_to_ascii_str(self.max_order, self._interval_set._intervals)
        elif format == 'json':
            return mocpy.spatial_moc_to_json_str(self.max_order, self._interval_set._intervals)
        else:
            formats = ('ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))

    @classmethod
    def from_string(cls, value, format='ascii'):
        """
        Deserialize the Spatial MOC from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.
        
        WARNING: the serialization must be strict, i.e. **must not** contain overlapping elements

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == 'ascii':
            intervals = mocpy.spatial_moc_from_ascii_str(value)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        elif format == 'json':
            intervals = mocpy.spatial_moc_from_json_str(value)
            return cls(IntervalSet(intervals, make_consistent=False), make_consistent=False)
        else:
            formats = ('ascii', 'json')
            raise ValueError('format should be one of %s' % (str(formats)))