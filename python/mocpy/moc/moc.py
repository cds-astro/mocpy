import contextlib
import functools
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from astropy import units as u
from astropy import wcs
from astropy.coordinates import (
    ICRS,
    Angle,
    BaseCoordinateFrame,
    Galactic,
    Latitude,
    Longitude,
    SkyCoord,
)
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file

with contextlib.suppress(ImportError):
    import regions
    from astropy_healpix import HEALPix

from .. import mocpy
from ..abstract_moc import AbstractMOC
from .boundaries import Boundaries
from .plot import border, fill
from .plot.wcs import WCS

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


def validate_lonlat(function):
    """Validate the longitude and latitudes entries of methods of the MOC class.

    Parameters
    ----------
    function : <class 'function'>
        must have the signature function(self, lon, lat, **kwargs)

    Returns
    -------
        applies desired transformations for the `lon` and `lat` arguments and calls
        `function` with these modified arguments

    Raises
    ------
    ValueError
        If `lon` and `lat` have inconsistent shapes.
    """

    @functools.wraps(function)
    def _validate_lonlat_wrap(self, lon, lat, **kwargs):
        # be sure that lon and lat are of the same shape
        if lon.shape != lat.shape:
            raise ValueError(
                f"'lon' and 'lat' should have the same shape but are of shapes {lon.shape} and {lat.shape}",
            )
        # lon and lat should be casted to arrays
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        # convert into astropy objects
        lon = lon if isinstance(lon, Longitude) else Longitude(lon)
        lat = lat if isinstance(lat, Latitude) else Latitude(lat)
        # convert to degrees
        lon = lon.to_value(u.deg).astype(np.float64)
        lat = lat.to_value(u.deg).astype(np.float64)
        return function(self, lon, lat, **kwargs)

    return _validate_lonlat_wrap


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
      python dictionary (see `from_json`).
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

    # Maximum order (or depth) of a MOC
    # (do not remove since it may be used externally).
    MAX_ORDER = np.uint8(29)

    def __init__(self, store_index):
        """Is a Spatial Coverage (S-MOC).

        Args:
            create_key: Object ensure __init__ is called by super-class/class-methods only
            store_index: index of the S-MOC in the rust-side storage
        """
        self.store_index = store_index

    @property
    def max_order(self):
        """Depth/order of the S-MOC."""
        depth = mocpy.get_smoc_depth(self.store_index)
        return np.uint8(depth)

    @classmethod
    def n_cells(cls, depth):
        """Get the number of cells for a given depth.

        Parameters
        ----------
        depth : int
            The depth. It is comprised between 0 and `~mocpy.moc.MOC.MAX_ORDER`

        Returns
        -------
        int
            The number of cells at the given order

        Examples
        --------
        >>> from mocpy import MOC
        >>> MOC.n_cells(0)
        12
        """
        if depth < 0 or depth > cls.MAX_ORDER:
            raise ValueError(
                f"The depth should be comprised between 0 and {cls.MAX_ORDER}, but {depth}"
                " was provided.",
            )
        return mocpy.n_cells_smoc(depth)

    def split_count(self, include_indirect_neighbours=False):
        """
        Return the number of disjoint MOCs the given MOC contains.

        Parameters
        ----------
        include_indirect_neighbours : bool
            if `false`, only consider  cells having a common edge as been part of a same MOC
            if `true`, also consider cells having a common vertex as been part of the same MOC
        """
        return mocpy.split_count(self.store_index, include_indirect_neighbours)

    def split(self, include_indirect_neighbours=False):
        """
        Return the disjoint MOCs this MOC contains.

        Parameters
        ----------
        include_indirect_neighbours : bool
            if `false`, only consider  cells having a common edge as been part of a same MOC
            if `true`, also consider cells having a common vertex as been part of the same MOC

        Warning
        -------
        Please use `~mocpy.moc.MOC.split_count` first to ensure the number is not too high
        """
        indices = mocpy.split(self.store_index, include_indirect_neighbours)
        return [MOC(index) for index in indices]

    def degrade_to_order(self, new_order):
        """
        Degrade the MOC instance to a new, less precise, MOC.

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
        index = mocpy.degrade(self.store_index, new_order)
        return MOC(index)

    def contains_skycoords(self, skycoords, keep_inside=True):
        """
        Return a boolean mask array of the positions lying inside (or outside) the MOC instance.

        Parameters
        ----------
        skycoords : `astropy.coordinates.SkyCoord`
            The sky coordinates that will be tested.
        keep_inside : bool, optional
            True by default. If so the mask describes coordinates lying inside the MOC. If ``keep_inside``
            is false, contains will return the mask of the coordinates lying outside the MOC.

        Returns
        -------
        array : `~np.ndarray`
            A mask boolean array

        See Also
        --------
        contains_lonlat
        """
        return self.contains_lonlat(
            lon=skycoords.icrs.ra,
            lat=skycoords.icrs.dec,
            keep_inside=keep_inside,
        )

    def contains(self, lon, lat, keep_inside=True):
        """Test wether a MOC contains --or not-- the given points. Returns a boolean mask array.

        .. deprecated:: 0.11.1
          `contains` is replaced by
          `contains_lonlat` for naming consistency.
          Please consider switching.

        Parameters
        ----------
        lon : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            Right ascension array in deg
        lat : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
            Declination array in deg
        keep_inside : bool, optional
            True by default. If so the mask describes coordinates lying inside the MOC. If ``keep_inside``
            is false, contains will return the mask of the coordinates lying outside the MOC.

        Returns
        -------
        array : `~np.ndarray`
            A mask boolean array

        See Also
        --------
        contains_skycoords
        """
        import warnings

        warnings.warn(
            "This method is deprecated and has been replaced by contains_lonlat",
            DeprecationWarning,
            stacklevel=2,
        )

        return self.contains_lonlat(lon, lat, keep_inside=keep_inside)

    @validate_lonlat
    def contains_lonlat(self, lon, lat, keep_inside=True):
        """Test wether a MOC contains (or not) the given points. Returns a boolean mask array.

        Parameters
        ----------
        lon : `astropy.coordinates.Longitude` or its supertype `astropy.units.Quantity`
            Right ascension array in deg
        lat : `astropy.coordinates.Latitude` or its supertype `astropy.units.Quantity`
            Declination array in deg
        keep_inside : bool, optional
            True by default. If so the mask describes coordinates lying inside the MOC. If ``keep_inside``
            is false, contains will return the mask of the coordinates lying outside the MOC.

        Raises
        ------
        ValueError : If `lon` and `lat` have mismatched shapes.

        Returns
        -------
        array : `~np.ndarray`
            A mask boolean array

        Examples
        --------
        >>> from mocpy import MOC
        >>> import numpy as np
        >>> from astropy.coordinates import Angle
        >>> import astropy.units as u
        >>> # create lists of coordinates
        >>> lon = Angle(np.array([[1, 2, 3], [-2, -40, -5]]), unit=u.deg)
        >>> lat = Angle(np.array([[20, 25, 10], [-60, 80, 0]]), unit=u.deg)
        >>> # create a polygonal moc from these
        >>> moc = MOC.from_polygon(lon=lon, lat=lat, max_depth=12)
        >>> moc.contains_lonlat(lon=lon, lat=lat) # returns all true
        array([[ True,  True,  True],
               [ True,  True,  True]])

        See Also
        --------
        contains_skycoords
        """
        mask = mocpy.filter_pos(
            self.store_index,
            lon,
            lat,
        )
        if keep_inside:
            return mask
        else:  # noqa: RET505
            return ~mask

    # TODO: implement: def contains_including_surrounding(self, lon, lat, distance)

    def fill(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Draw the MOC on a matplotlib axis.

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
        >>> from mocpy import MOC
        >>> import astropy.units as u
        >>> import matplotlib.pyplot as plt
        >>> # Create a MOC
        >>> moc = MOC.from_ring(external_radius=2*u.deg,
        ...                     internal_radius=1*u.deg,
        ...                     lat=0*u.rad, lon=0*u.rad,
        ...                     max_depth=13,
        ...                    )
        >>> # Plot the MOC using matplotlib
        >>> fig = plt.figure(figsize=(10, 10))
        >>> wcs = moc.wcs(fig)
        >>> ax = fig.add_subplot(projection=wcs)
        >>> moc.fill(ax, wcs, color='blue')
        """
        fill.fill(self, ax, wcs, **kw_mpl_pathpatch)

    def border(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Draw the MOC border.s on a matplotlib axis.

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
        >>> from mocpy import MOC
        >>> from astropy.coordinates import Latitude, Longitude
        >>> import astropy.units as u
        >>> import matplotlib.pyplot as plt
        >>> # Create a MOC
        >>> lon = Longitude([5, -5, -5, 5], u.deg)
        >>> lat = Latitude([5, 5, -5, -5], u.deg)
        >>> moc = MOC.from_polygon(lon, lat)
        >>> # Plot the MOC using matplotlib
        >>> fig = plt.figure(figsize=(10, 10))
        >>> wcs = moc.wcs(fig)
        >>> ax = fig.add_subplot(projection=wcs)
        >>> moc.border(ax, wcs, color='blue')
        """
        border.border(self, ax, wcs, **kw_mpl_pathpatch)

    def get_boundaries(self, order=None):
        """
        Return the sky coordinates defining the border(s) of the MOC.

        The border(s) are expressed as a list of SkyCoord.
        Each SkyCoord refers to the coordinates of one border of the MOC (i.e.
        either a border of a connected MOC part or a border of a hole
        located in a connected MOC part).
        This function is currently not stable: encoding a vertex of a
        HEALPix cell (N, E, S, W) should not depend on the position of the
        vertex but rather on the uniq value (+ 2 bits to encode the direction
        of the vertex).

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

        warnings.warn(
            "This method is not stable. A future more stable algorithm will be implemented!",
            DeprecationWarning,
            stacklevel=2,
        )
        return Boundaries.get(self, order)

    @classmethod
    def from_fits_image(cls, hdu, max_norder, mask=None):
        """
        Create a `~mocpy.moc.MOC` from an image stored as a FITS file.

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

        # Compute a WCS from the header of the image
        w = wcs.WCS(header)

        if mask is None:
            data = hdu.data
            # A mask is computed discarding nan floating values
            mask = np.isfinite(data)

            # If the BLANK keyword is set to a value then we mask those
            # pixels too
            if header.get("BLANK") is not None:
                discard_val = header["BLANK"]

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
        skycrd = SkyCoord(world, unit="deg", frame=frame)

        # Compute the order based on the CDELT
        # Should work even if the fits is defined using cd_ij or crota, see:
        # https://docs.astropy.org/en/stable/api/astropy.wcs.Wcsprm.html#astropy.wcs.Wcsprm.get_cdelt
        [c1, c2] = w.wcs.get_cdelt()

        max_res_px = np.sqrt(c1 * c1 + c2 * c2) * np.pi / 180.0
        max_depth_px = int(np.floor(np.log2(np.pi / (3 * max_res_px * max_res_px)) / 2))

        max_norder = min(max_norder, max_depth_px)

        moc = MOC.from_lonlat(
            lon=skycrd.icrs.ra,
            lat=skycrd.icrs.dec,
            max_norder=max_norder,
        )
        return moc  # noqa: RET504

    @classmethod
    def from_fits_images(cls, path_l, max_norder, hdu_index=0):
        """
        Load a MOC from a set of FITS file images.

        Assumes the data of the image is stored in the first HDU of the FITS file.
        Please call `~mocpy.moc.MOC.from_fits_image` for passing another hdu than the first one.

        Parameters
        ----------
        path_l : [str]
            A list of path where the fits image are located.
        max_norder : int
            The MOC resolution.
        hdu_index : int
            Index of the the HDU containing the image in each FITS file (default = 0)

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The union of all the MOCs created from the paths found in ``path_l``.
        """
        moc = MOC.new_empty(max_norder)
        for filename in path_l:
            with fits.open(filename) as hdul:
                current_moc = MOC.from_fits_image(
                    hdu=hdul[hdu_index],
                    max_norder=max_norder,
                )
                moc = moc.union(current_moc)

        return moc

    @classmethod
    def from_vizier_table(cls, table_id, nside=256):
        """
        Create a `~mocpy.moc.MOC` object from a VizieR table.

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
            raise ValueError(
                f"Bad value for nside. Must be in {nside_possible_values}",
            )
        return cls.from_ivorn("ivo://CDS/" + table_id, nside)

    MOC_SERVER_ROOT_URL = "http://alasky.unistra.fr/MocServer/query"

    @classmethod
    def from_ivorn(cls, ivorn, nside=256):
        """
        Create a `~mocpy.moc.MOC` object from a given ivorn.

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
        return cls.from_url(
            "{}?{}".format(
                MOC.MOC_SERVER_ROOT_URL,
                urlencode({"ivorn": ivorn, "get": "moc", "order": int(np.log2(nside))}),
            ),
        )

    @classmethod
    def from_url(cls, url):
        """
        Create a `~mocpy.moc.MOC` object from a given url.

        Parameters
        ----------
        url : str
            The url of a FITS file storing a MOC.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC.
        """
        # TODO: as is, this is a duplicate of abstract class `from_fits` called with an url
        path = download_file(url, show_progress=False, timeout=60)
        return cls.load(path, "fits")

    @classmethod
    def from_skycoords(cls, skycoords, max_norder):
        """
        Create a MOC from an `astropy.coordinates.SkyCoord`.

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
        return cls.from_lonlat(
            lon=skycoords.icrs.ra,
            lat=skycoords.icrs.dec,
            max_norder=max_norder,
        )

    @classmethod
    @validate_lonlat
    def from_lonlat(cls, lon, lat, max_norder):
        """
        Create a MOC from astropy lon, lat `astropy.units.Quantity`.

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
        index = mocpy.from_lonlat(
            max_norder,
            lon,
            lat,
        )
        return cls(index)

    @classmethod
    def from_multiordermap_fits_file(
        cls,
        path,
        cumul_from=0.0,
        cumul_to=1.0,
        asc=False,
        strict=True,
        no_split=True,
        reverse_decent=False,
    ):
        """
        Create a MOC from a mutli-order map FITS file.

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
        path : str or pathlib.Path
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
        index = mocpy.spatial_moc_from_multiordermap_fits_file(
            str(path),
            np.float64(cumul_from),
            np.float64(cumul_to),
            asc,
            strict,
            no_split,
            reverse_decent,
        )
        return cls(index)

    @classmethod
    def probabilities_in_multiordermap(cls, mocs, multiordermap, n_threads=None):
        """Calculate the probabilities in the intersection between the multiordermap and the MOCs.

        Multi-MOC version of `probability_in_multiordermap`. This is parallelized.

        Parameters
        ----------
        mocs : list[mocpy.MOC]
            A list of `mocpy.MOC`.
        multiordermap : astropy.table.Table, or astropy.table.QTable
            Should have a column ``UNIQ`` that
            corresponds to HEALPix cells and a ``PROBDENSITY`` column.
        n_threads : int
            Number of threads to be used (all available threads by default)

        Returns
        -------
        list[float]
            A list containing the probability for each MOC

        Notes
        -----
        In wasm compilations (ex for pyodide), this won't raise an error, but will be
        single-threaded.

        """
        if not isinstance(multiordermap, Table):
            raise ValueError(
                "An argument of type 'astropy.table.Table'"
                f" is expected. Got '{type(multiordermap)}'",
            )

        indices = np.array(
            [moc.store_index for moc in mocs],
            dtype=cls._store_index_dtype(),
        )

        uniq = multiordermap["UNIQ"]
        probdensity = multiordermap["PROBDENSITY"]

        try:
            uniq_mask = uniq.data.mask
            probdensity_mask = probdensity.data.mask
        except AttributeError:
            uniq_mask = np.zeros(uniq.shape)
            probdensity_mask = np.zeros(probdensity.shape)

        return mocpy.multi_multiorder_probdens_map_sum_in_smoc(
            indices,
            np.array(uniq.data, dtype="uint64"),
            np.array(uniq_mask, dtype="bool"),
            np.array(probdensity.data, dtype="float"),
            np.array(probdensity_mask, dtype="bool"),
            n_threads,
        )

    def probability_in_multiordermap(self, multiordermap):
        """Calculate the probability in the intersection between the multiordermap and the MOC.

        ``PROBDENSITY`` values are multiplied by the area of their associated HEALPix
        cell before summing them. For cells that are not complete, the ratio of the area
        is used.

        Parameters
        ----------
        multiordermap : str, pathlib.Path, astropy.table.Table, or astropy.table.QTable
            If ``multiordermap`` is given as a string or `~pathlib.Path`, the probability
            will be read from the column ``PROBDENSITY`` of the FITS file.
            If it is an `~astropy.table.Table`, then it should have a column ``UNIQ`` that
            corresponds to HEALPix cells and a ``PROBDENSITY`` column.

        Returns
        -------
        float
            The probability in the intersection between the MOC and the Multi-Order-Map
            coverages.

        Examples
        --------
        >>> from mocpy import MOC
        >>> import numpy as np
        >>> from astropy.table import Table
        >>> all_sky = MOC.from_str("0/0-11")
        >>> # Let's create a meaningless multiorder map
        >>> uniq = [4 * 4**4 + x for x in range(20)]
        >>> rng = np.random.default_rng(0)
        >>> probdensity = rng.random(20) / 100
        >>> multi_order_map = Table([uniq, probdensity], names=("UNIQ", "PROBDENSITY"))
        >>> # The probability to be in the intersection with the all sky is
        >>> round(all_sky.probability_in_multiordermap(multi_order_map), 4)
        0.0004

        See Also
        --------
        probabilities_in_multiordermap: makes this calculation for a list of MOCs in a
        parallelized way.

        """
        index = self.store_index

        if isinstance(multiordermap, Table):
            uniq = multiordermap["UNIQ"]
            probdensity = multiordermap["PROBDENSITY"]
            try:
                uniq_mask = uniq.data.mask
                probdensity_mask = probdensity.data.mask
            except AttributeError:
                uniq_mask = np.zeros(uniq.shape)
                probdensity_mask = np.zeros(probdensity.shape)

            return mocpy.multiorder_probdens_map_sum_in_smoc(
                index,
                np.array(uniq.data, dtype="uint64"),
                np.array(uniq_mask, dtype="bool"),
                np.array(probdensity.data, dtype="float"),
                np.array(probdensity_mask, dtype="bool"),
            )
        if isinstance(multiordermap, (Path, str)):
            return mocpy.multiordermap_sum_in_smoc_from_file(index, str(multiordermap))

        raise ValueError(
            "An argument of type 'str', 'pathlib.Path', or 'astropy.table.Table'"
            f" is expected. Got '{type(multiordermap)}'",
        )

    def sum_in_multiordermap(self, multiordermap: Table, column: str):
        """Calculate the sum of a column from a multiordermap in the intersection with the MOC.

        Parameters
        ----------
        multiordermap : astropy.table.Table
            The table should have a ``UNIQ`` that corresponds to HEALPix cells in the uniq
            notation.
        column : str
            The name of the column to sum.

        Returns
        -------
        float
            The sum of the values in the intersection between the MOC and the
            multiorder map coverages.

        Examples
        --------
        >>> from mocpy import MOC
        >>> import numpy as np
        >>> from astropy.table import Table
        >>> all_sky = MOC.from_str("0/0-11")
        >>> # Let's create a meaningless multiorder map
        >>> uniq = [4 * 4**5 + x for x in np.arange(200)]
        >>> rng = np.random.default_rng(0)
        >>> column = rng.random(200)
        >>> multi_order_map = Table([uniq, column], names=("UNIQ", "column"))
        >>> round(all_sky.sum_in_multiordermap(multi_order_map, "column"), 4)
        107.9259

        """
        index = self.store_index
        uniq = multiordermap["UNIQ"]
        values_to_sum = multiordermap[column]
        try:
            uniq_mask = uniq.data.mask
            values_to_sum_mask = values_to_sum.data.mask
        except AttributeError:
            uniq_mask = np.zeros(uniq.shape)
            values_to_sum_mask = uniq_mask

        return mocpy.multiordermap_sum_in_smoc(
            index,
            np.array(uniq.data, dtype="uint64"),
            np.array(uniq_mask, dtype="bool"),
            np.array(values_to_sum.data, dtype="float"),
            np.array(values_to_sum_mask, dtype="bool"),
        )

    @classmethod
    def from_valued_healpix_cells(
        cls,
        uniq,
        values,
        max_depth=None,
        values_are_densities=False,
        cumul_from=0.0,
        cumul_to=1.0,
        asc=False,
        strict=True,
        no_split=True,
        reverse_decent=False,
    ):
        """
        Create a MOC from a list of uniq associated with values.

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
            Value associated with each ``uniq`` cells. dtype must be np.float64
        max_depth : int
            The max depth of the MOC, should be at least as large as the depth corresponding of the smallest HEALPix cell found in ``uniq``.
            Warnings:
            1 - the depth of the returned MOC will be at least as deep as the smallest HEALPix cell found in ``uniq``.
            2 - contrary to MOCPy before v0.12, the user has to degrade the MOC if `max_depth` < smallest HEALPix cell depth.
        values_are_densities: tell whether the values depend on the cell area or not
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
        if max_depth is None:
            import warnings

            warnings.warn(
                "To avoid an extra loop, it is preferable to provide the max_depth parameter."
                "It will probably become mandatory in future releases.",
                UserWarning,
                stacklevel=2,
            )

            max_depth = int(np.log2(uniq.max() >> 2)) >> 1
            if max_depth < 0 or max_depth > 29:
                raise ValueError(
                    "Invalid uniq numbers. Too big uniq or negative uniq numbers might be the cause.",
                )

        index = mocpy.from_valued_hpx_cells(
            np.uint8(max_depth),
            uniq.astype(np.uint64),
            values.astype(np.float64),
            values_are_densities,
            np.float64(cumul_from),
            np.float64(cumul_to),
            asc,
            strict,
            no_split,
            reverse_decent,
        )
        return cls(index)

    @classmethod
    @validate_lonlat
    def from_elliptical_cone(cls, lon, lat, a, b, pa, max_depth, delta_depth=2):
        """
        Create a MOC from an elliptical cone.

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
        index = mocpy.from_elliptical_cone(
            lon[0],
            lat[0],
            np.float64(a.to_value(u.deg)),
            np.float64(b.to_value(u.deg)),
            np.float64(pa.to_value(u.deg)),
            np.uint8(max_depth),
            np.uint8(delta_depth),
        )
        return cls(index)

    @classmethod
    @validate_lonlat
    def from_cone(cls, lon, lat, radius, max_depth, delta_depth=2):
        """
        Create a MOC from a cone.

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
        index = mocpy.from_cone(
            lon[0],
            lat[0],
            np.float64(radius.to_value(u.deg)),
            np.uint8(max_depth),
            np.uint8(delta_depth),
        )
        return cls(index)

    @classmethod
    @validate_lonlat
    def from_box(cls, lon, lat, a, b, angle, max_depth):
        """
        Create a MOC from a box/rectangle.

        The box is centered around the (`lon`, `lat`) position.

        Parameters
        ----------
        lon : `~astropy.coordinates.Longitude` or its supertype `~astropy.units.Quantity`
            The longitude of the center of the cone.
        lat : `~astropy.coordinates.Latitude` or its supertype `~astropy.units.Quantity`
            The latitude of the center of the cone.
        a : `~astropy.coordinates.Angle`
            The semi-major axis of the box, i.e. half of the box's width.
        b : `~astropy.coordinates.Angle`
            The semi-minor axis of the box, i.e. half of the box's height.
        angle : `astropy.coordinates.Angle`
            The tilt angle between the north and the semi-major axis, east of north.
        max_depth : int
            Maximum HEALPix cell resolution.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC

        Examples
        --------
        >>> from mocpy import MOC
        >>> from astropy.coordinates import Angle, Longitude, Latitude
        >>> moc = MOC.from_box(
        ...  lon=Longitude("0d"),
        ...  lat=Latitude("0d"),
        ...  a=Angle("10d"),
        ...  b=Angle("2d"),
        ...  angle=Angle("30d"),
        ...  max_depth=10
        ... )
        """
        index = mocpy.from_box(
            lon[0],
            lat[0],
            np.float64(a.to_value(u.deg)),
            np.float64(b.to_value(u.deg)),
            np.float64(angle.to_value(u.deg)),
            np.uint8(max_depth),
        )
        return cls(index)

    @classmethod
    @validate_lonlat
    def from_ring(
        cls,
        lon,
        lat,
        internal_radius,
        external_radius,
        max_depth,
        delta_depth=2,
    ):
        """
        Create a MOC from a ring.

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
        index = mocpy.from_ring(
            lon[0],
            lat[0],
            np.float64(internal_radius.to_value(u.deg)),
            np.float64(external_radius.to_value(u.deg)),
            np.uint8(max_depth),
            np.uint8(delta_depth),
        )
        return cls(index)

    @classmethod
    def from_polygon_skycoord(cls, skycoord, max_depth=10):
        """
        Create a MOC from a polygon.

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
        return cls.from_polygon(
            lon=skycoord.icrs.ra,
            lat=skycoord.icrs.dec,
            max_depth=np.uint8(max_depth),
        )

    @classmethod
    def from_polygons(cls, list_vertices, max_depth=10, n_threads=None):
        """
        Create a list of MOCs list from a list of polygons.

        Parameters
        ----------
        list_vertices : list[`~astropy.coordinates.SkyCoord`]
            A list of `~astropy.coordinates.SkyCoord` objects. Each SkyCoord object should
            contain more than three vertices and they should each describe a polygon.
         max_depth : int, optional
            The resolution of the MOC. Set to 10 by default.
        n_threads : int, optional
            The number of threads to be used

        Returns
        -------
        list[`mocpy.MOC`]

        Examples
        --------
        >>> from astropy.coordinates import SkyCoord
        >>> from mocpy import MOC
        >>> list_vertices = [
        ...     SkyCoord([-4, 4, 4, -4], [4, 4, -4, -4], unit="deg"),
        ...     SkyCoord([0, 6, 0, -6], [6, 0, -6, 0], unit="deg")
        ... ]
        >>> list_mocs = MOC.from_polygons(list_vertices)
        >>> len(list_mocs) # two mocs are created
        2

        """
        lon_lat_list = [
            f(x)
            for x in list_vertices
            for f in (lambda x: x.icrs.ra.deg, lambda x: x.icrs.dec.deg)
        ]
        indices = mocpy.from_polygons(lon_lat_list, np.uint8(max_depth), n_threads)
        return [cls(index) for index in indices]

    @classmethod
    @validate_lonlat
    def from_polygon(cls, lon, lat, complement=False, max_depth=10):
        """
        Create a MOC from a polygon.

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
        complement : return the complement of the polygon. Set to False by default.
        max_depth : int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The resulting MOC
        """
        index = mocpy.from_polygon(
            lon,
            lat,
            complement,
            np.uint8(max_depth),
        )
        return cls(index)

    @classmethod
    def from_stcs(cls, stcs, max_depth, delta_depth=2):
        """
        Create a MOC from a STC-S.

        Info
        ----
        Time, Spectral and Redshift sub-phrases are ignored.

        Warning
        -------
        There is so far no implicit conversion, so the STC-S string will be rejected if:
        * the frame is different from `ICRS`
        * the flavor is different from `Spher2`
        * the units are different from `degrees`
        The implementation is not (yet?) fully compliant with the STC standard (see MOC Lib rust for more details):
        * DIFFERENCE is so far interpreted as a symmetrical difference (XOR) while it is a MINUS in the STC standard
        * Self-intersecting Polygons are supported, and the "interior" usually has the smallest area

        Parameters
        ----------
        stcs : str
            The STC-S string.
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
        >>> moc1 = MOC.from_stcs("Circle ICRS 147.6 69.9 0.4", max_depth=14)
        >>> moc2 = MOC.from_cone(lon=147.6 * u.deg, lat=69.9 * u.deg, radius=Angle(0.4, u.deg), max_depth=14)
        >>> assert (moc1 == moc2)
        """
        index = mocpy.from_stcs(stcs, np.uint8(max_depth), np.uint8(delta_depth))
        return cls(index)

    @classmethod
    def from_astropy_regions(cls, region, max_depth: int):
        """Create a SMOC from an astropy regions.

        This creates the MOC of the requested order that contains entirely the astropy
        region. See https://github.com/astropy/regions.

        Parameters
        ----------
        region : `~regions.SkyRegion`
            The supported sky regions are `~regions.CircleSkyRegion`,
            `~regions.CircleAnnulusSkyRegion`, `~regions.EllipseSkyRegion`,
            `~regions.RectangleSkyRegion`, `~regions.PolygonSkyRegion`,
            `~regions.PointSkyRegion`.
        max_depth : int
            The maximum HEALPix cell resolution of the MOC. Should be comprised between
            0 and 29.

        Returns
        -------
        `~mocpy.moc.MOC`

        Examples
        --------
        >>> from astropy.coordinates import SkyCoord
        >>> point = SkyCoord("+23h20m48.3s +61d12m06s")
        >>> point_region = regions.PointSkyRegion(point)
        >>> moc_point = MOC.from_astropy_regions(point_region, max_depth=10)
        >>> moc_point
        10/3663728
        """
        supported_regions_types = (
            regions.CircleSkyRegion,
            regions.CircleAnnulusSkyRegion,
            regions.EllipseSkyRegion,
            regions.RectangleSkyRegion,
            regions.PolygonSkyRegion,
            regions.PointSkyRegion,
        )
        if isinstance(region, regions.CircleSkyRegion):
            center = region.center.icrs
            return cls.from_cone(
                center.ra,
                center.dec,
                radius=region.radius,
                max_depth=max_depth,
            )
        if isinstance(region, regions.CircleAnnulusSkyRegion):
            center = region.center.icrs
            return cls.from_ring(
                center.ra,
                center.dec,
                internal_radius=region.inner_radius,
                external_radius=region.outer_radius,
                max_depth=max_depth,
            )
        if isinstance(region, regions.EllipseSkyRegion):
            center = region.center.icrs
            return cls.from_elliptical_cone(
                center.ra,
                center.dec,
                a=region.height / 2,
                b=region.width / 2,
                pa=region.angle,
                max_depth=max_depth,
            )
        if isinstance(region, regions.RectangleSkyRegion):
            center = region.center.icrs
            return cls.from_box(
                center.ra,
                center.dec,
                a=region.width / 2,
                b=region.height / 2,
                angle=region.angle + Angle("90d"),
                max_depth=max_depth,
            )
        if isinstance(region, regions.PolygonSkyRegion):
            return cls.from_polygon_skycoord(region.vertices, max_depth=max_depth)
        if isinstance(region, regions.PointSkyRegion):
            return cls.from_skycoords(region.center, max_norder=max_depth)

        raise ValueError(
            "'from_astropy_regions' does not support this region type."
            f"The supported regions are: {supported_regions_types}",
        )

    @classmethod
    def new_empty(cls, max_depth):
        """
        Create a new empty MOC of given depth.

        Parameters
        ----------
        max_depth : int, The resolution of the MOC


        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The MOC
        """
        index = mocpy.new_empty_smoc(np.uint8(max_depth))
        return cls(index)

    @classmethod
    def from_healpix_cells(cls, ipix, depth, max_depth):
        """
        Create a MOC from a set of HEALPix cells at various depths.

        Parameters
        ----------
        ipix : `numpy.ndarray`
            HEALPix cell indices in the NESTED notation. dtype must be np.uint64
        depth : `numpy.ndarray`
            Depth of the HEALPix cells. Must be of the same size of `ipix`.
            dtype must be np.uint8. Corresponds to the `level` of an HEALPix cell in astropy.healpix.
        max_depth : int, The resolution of the MOC (degrades on the fly input cells if necessary)

        Raises
        ------
        IndexError
            When `ipix` and `depth` do not have the same shape

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The MOC
        """
        index = mocpy.from_healpix_cells(
            np.uint8(max_depth),
            depth.astype(np.uint8),
            ipix.astype(np.uint64),
        )
        return cls(index)

    @classmethod
    def from_depth29_ranges(cls, max_depth, ranges):
        """
        Create a MOC from a set of ranges of HEALPix Nested indices at order 29.

        For each range, the lower bound is inclusive and the upper bound is exclusive.
        The range `[0, 12*4^29[` represents the full-sky.

        Parameters
        ----------
        max_depth : int, The resolution of the MOC
        ranges : `~numpy.ndarray`, optional
                a N x 2 numpy array representing the set of depth 29 HEALPix nested ranges.
                defaults to `np.zeros((0, 2), dtype=np.uint64)`

        Returns
        -------
        moc : `~mocpy.moc.MOC`
            The MOC
        """
        ranges = np.zeros((0, 2), dtype=np.uint64) if ranges is None else ranges

        if ranges.shape[1] != 2:
            raise ValueError("expected a N x 2 numpy array for ranges")

        if ranges.dtype is not np.uint64:
            ranges = ranges.astype(np.uint64)

        index = mocpy.from_hpx_ranges(np.uint8(max_depth), ranges)
        return cls(index)

    @classmethod
    def from_stmoc_time_fold(cls, tmoc, stmoc):
        """
        Build a new S-MOC from the fold operation of the given ST-MOC by the given T-MOC.

        Parameters
        ----------
        tmoc : `~mocpy.tmoc.TimeMoc`
        stmoc : `~mocpy.stmoc.STMoc`
        """
        store_index = mocpy.project_on_second_dim(tmoc.store_index, stmoc.store_index)
        return cls(store_index)

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
        return Angle(np.sqrt(np.pi / (3 * 4 ** (order))), unit="rad")

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
        order = np.ceil(np.log2(np.pi / (3 * res_rad * res_rad)) / 2)
        return np.uint8(order)

    @property
    def _fits_header_keywords(self):
        return {
            "PIXTYPE": "HEALPIX",
            "ORDERING": "NUNIQ",
            "COORDSYS": ("C", "reference frame (C=ICRS)"),
            "MOCORDER": self.max_order,
            "MOCTOOL": "MOCPy",
        }

    @property
    def _fits_format(self):
        depth = self.max_order
        return "1J" if depth <= 13 else "1K"

    @property
    def sky_fraction(self):
        """Sky fraction covered by the MOC.

        Examples
        --------
        >>> from mocpy import MOC
        >>> MOC.from_string("0/0-11").sky_fraction # all sky
        1.0
        """
        return mocpy.coverage_fraction(self.store_index)

    # TODO : move this in astroquery.Simbad.query_region
    # See https://github.com/astropy/astroquery/pull/1466
    def query_simbad(self, max_rows=10000, timeout=1000):
        """Query a view of SIMBAD data for SIMBAD objects in the coverage of the MOC instance.

        Parameters
        ----------
        max_rows : int, optional
                maximum number of row returned
        timeout : float, optional
                timeout before aborting the query, default to 1000s
        """
        return self._query("SIMBAD", max_rows, timeout)

    # TODO : move this in astroquery.Vizier.query_region
    # See https://github.com/astropy/astroquery/pull/1466
    def query_vizier_table(self, table_id: str, max_rows=10000, timeout=1000):
        """Query a VizieR table for sources in the coverage of the MOC instance.

        Parameters
        ----------
        table_id : str
                corresponds to a VizieR table id
        max_rows : int, optional
                maximum number of row returned
        timeout : float, optional
                timeout before aborting the query, default to 1000s
        """
        return self._query(table_id, max_rows, timeout)

    # TODO : move this in astroquery
    def _query(self, resource_id, max_rows=100000, timeout=1000):
        """
        Query Simbad or a VizieR table.

        Find sources in the coverage of the MOC instance.
        """
        import requests
        from astropy.io.votable import parse_single_table

        moc_file = BytesIO()
        moc_fits = self.serialize(format="fits", pre_v2=True)
        moc_fits.writeto(moc_file)

        r = requests.post(
            "http://cdsxmatch.u-strasbg.fr/QueryCat/QueryCat",
            data={
                "mode": "mocfile",
                "catName": resource_id,
                "format": "votable",
                "limit": max_rows,
            },
            files={"moc": moc_file.getvalue()},
            headers={"User-Agent": "MOCPy"},
            stream=True,
            timeout=timeout,
        )

        votable = BytesIO()
        votable.write(r.content)

        return parse_single_table(votable).to_table()

    def wcs(
        self,
        fig,
        coordsys="icrs",
        projection="AIT",
        rotation=None,
    ):
        """
        Get a wcs that can be given to matplotlib to plot the MOC.

        Parameters
        ----------
        fig : `~matplotlib.pyplot.figure`
            The matplotlib figure used for plotting the MOC.
        coordsys : str, optional
            Coordinate system. Default to "icrs". Must be in ["icrs", "galactic"].
        projection : str, optional
            World base -> Image base projection type.
            See http://docs.astropy.org/en/stable/wcs/#supported-projections for
            the projections currently supported in astropy. Default to Aitoff.
        rotation : `~astropy.coordinates.Angle`, optional
            The angle of rotation. Default to no rotation.

        Returns
        -------
        wcs : `~astropy.wcs.WCS`
            The WCS that can be passed to mocpy.MOC.fill/border.

        Examples
        --------
        >>> from mocpy import MOC
        >>> import matplotlib.pyplot as plt
        >>> moc = MOC.from_str("2/2-25 28 29 4/0 6/")
        >>> fig = plt.figure()
        >>> moc.wcs(fig)
        WCS Keywords
        <BLANKLINE>
        Number of WCS axes: 2
        CTYPE : 'RA---AIT'  'DEC--AIT'
        CRVAL : 92.29966711743452  54.33295312309193
        CRPIX : 320.0  240.0
        PC1_1 PC1_2  : 1.0  -0.0
        PC2_1 PC2_2  : 0.0  1.0
        CDELT : -0.27794934051515896  0.27794934051515896
        NAXIS : 0  0
        """
        if rotation is None:
            rotation = Angle(0, u.radian)
        # The center is set to the barycenter of all its HEALPix cells
        center = self.barycenter()
        # The fov is computed from the largest distance between the center and any cells of it
        fov = 2 * self.largest_distance_from_coo_to_vertices(center)

        return WCS(
            fig,
            fov=fov,
            center=center,
            coordsys=coordsys,
            rotation=rotation,
            projection=projection,
        ).w

    def plot(self, title="MOC", frame=None):
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

        warnings.warn(
            "This method is deprecated and is no longer tested."
            "Please refer to `MOC.fill` and `MOC.border` methods",
            DeprecationWarning,
            stacklevel=2,
        )

        frame = ICRS() if frame is None else frame

        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap

        plot_order = 8
        plotted_moc = (
            self.degrade_to_order(plot_order) if self.max_order > plot_order else self
        )

        num_pixels_map = 1024
        delta = 2.0 * np.pi / num_pixels_map

        x = np.arange(-np.pi, np.pi, delta)
        y = np.arange(-np.pi / 2, np.pi / 2, delta)
        lon_rad, lat_rad = np.meshgrid(x, y)
        hp = HEALPix(nside=(1 << plotted_moc.max_order), order="nested")

        if frame and not isinstance(frame, BaseCoordinateFrame):
            raise ValueError(
                "Only Galactic/ICRS coordinate systems are supported."
                "Please set `coord` to either 'C' or 'G'.",
            )

        pix_map = hp.lonlat_to_healpix(lon_rad * u.rad, lat_rad * u.rad)

        m = np.zeros(12 * 4 ** (plotted_moc.max_order))
        pix_id = plotted_moc.flatten()

        # change the HEALPix cells if the frame of the MOC is not the same as the one associated with the plot method.
        if isinstance(frame, Galactic):
            lon, lat = hp.boundaries_lonlat(pix_id, step=2)
            sky_crd = SkyCoord(lon, lat, unit="deg")
            pix_id = hp.lonlat_to_healpix(sky_crd.galactic.l, sky_crd.galactic.b)

        m[pix_id] = 1

        z = np.flip(m[pix_map], axis=1)

        plt.figure(figsize=(10, 10))

        ax = plt.subplot(111, projection="mollweide")
        ax.set_xticklabels(
            [
                "150°",
                "120°",
                "90°",
                "60°",
                "30°",
                "0°",
                "330°",
                "300°",
                "270°",
                "240°",
                "210°",
                "180°",
            ],
        )

        color_map = LinearSegmentedColormap.from_list("w2r", ["#eeeeee", "#aa0000"])
        color_map.set_under("w")
        color_map.set_bad("gray")

        ax.pcolormesh(x, y, z, cmap=color_map, vmin=0, vmax=1)
        ax.tick_params(labelsize=14, labelcolor="#000000")
        plt.title(title)
        plt.grid(visible=True, linestyle="--", linewidth=1, color="#555555")

        plt.show()

    @classmethod
    def load(cls, path, format="fits"):  # noqa: A002
        """
        Load the Spatial MOC from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str or pathlib.Path
            The path to the file to load the MOC from.
        format : str, optional
            The format from which the MOC is loaded.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        """
        path = str(path)
        if format == "fits":
            index = mocpy.spatial_moc_from_fits_file(path)
            return cls(index)
        if format == "ascii":
            index = mocpy.spatial_moc_from_ascii_file(path)
            return cls(index)
        if format == "json":
            index = mocpy.spatial_moc_from_json_file(path)
            return cls(index)
        formats = ("fits", "ascii", "json")
        raise ValueError("format should be one of %s" % (str(formats)))

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load MOC from raw bytes of a FITS file."""
        index = mocpy.spatial_moc_from_fits_raw_bytes(raw_bytes)
        return cls(index)

    @classmethod
    def from_string(cls, value, format="ascii"):  # noqa: A002
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
        if format == "ascii":
            index = mocpy.spatial_moc_from_ascii_str(value)
            return cls(index)
        if format == "json":
            index = mocpy.spatial_moc_from_json_str(value)
            return cls(index)
        formats = ("ascii", "json")
        raise ValueError("format should be one of %s" % (str(formats)))

    @property
    def uniq_hpx(self):
        """
        Return a `np.array` of the uniq HEALPIx indices of the cell in the MOC.

        Warning
        -------
        The output is not sorted, the order follow the order of HEALPix cells in
        the underlying sorted array of depth29 nested ranges, i.e. the natural order
        of the cells is the underlying z-order curve.
        """
        return mocpy.to_uniq_hpx(self.store_index)

    @property
    def to_depth29_ranges(self):
        """Return the list of order 29 HEALPix nested ranges this MOC contains."""
        return mocpy.to_ranges(self.store_index)

    def to_rgba(self, y_size=300):
        """
        Create a matplotlib compatible RGBA preview of the given MOC.

        Parameters
        ----------
        y_size : the number of pixels along the y-axis

        Returns
        -------
        array : A (2 * y_size, y_size, 4) array of 0-255 int values.
        """
        return mocpy.to_rgba(self.store_index, y_size)

    def display_preview(self, y_size=300):
        """
        Display a preview of the MOC (calling internally the `to_rgba` method).

        Parameters
        ----------
        y_size : the number of pixels along the y-axis, default value is 300
        """
        import matplotlib.pyplot as plt

        plt.axis("off")
        plt.imshow(mocpy.to_rgba(self.store_index, y_size))
        plt.show()

    def barycenter(self):
        """Return the Barycenter of the MOC."""
        lonlat = mocpy.get_barycenter(self.store_index)
        return SkyCoord(lonlat[0], lonlat[1], unit="rad")

    def largest_distance_from_coo_to_vertices(self, coo):
        """Return the largest distance between the given coordinates and vertices of the MOC cells."""
        angle = mocpy.get_largest_distance_from_coo_to_moc_vertices(
            self.store_index,
            coo.ra.rad,
            coo.dec.rad,
        )
        return angle * u.rad
