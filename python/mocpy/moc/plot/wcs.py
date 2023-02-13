import numpy as np

from astropy import coordinates
from astropy import wcs

import astropy.units as u


class WCS:
    """
    Create a WCS for vizualizing a MOC in a matplotlib axis.

    Parameters
    ----------
    fig : `~matplotlib.pyplot.figure`
        The matplotlib figure used for plotting the MOC.
    fov : `~astropy.units.Quantity`
        Size of the field of view.
    center : `~astropy.coordinates.SkyCoord`, optional
        World coordinates matching with the center of the plot. Default to (0 deg, 0 deg) (in ICRS frame).
    coordsys : str, optional
        Coordinate system. Default to "icrs". Must be in ["icrs", "galactic"].
    projection : str, optional
        World base -> Image base projection type. See http://docs.astropy.org/en/stable/wcs/#supported-projections for
        the projections currently supported in astropy. Default to Aitoff.
    rotation : `~astropy.coordinates.Angle`, optional
        The angle of rotation. Default to no rotation.

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        The WCS that can be passed to mocpy.MOC.fill/border.

    Examples
    --------
    >>> from mocpy import MOC, WCS
    >>> from astropy.coordinates import Angle, SkyCoord
    >>> import astropy.units as u
    >>> # Load a MOC
    >>> filename = './../resources/P-GALEXGR6-AIS-FUV.fits'
    >>> moc = MOC.from_fits(filename)
    >>> # Plot the MOC using matplotlib
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure(111, figsize=(15, 15))
    >>> # Define a WCS as a context
    >>> with WCS(fig,
    ...         fov=200 * u.deg,
    ...         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
    ...         coordsys="icrs",
    ...         rotation=Angle(0, u.degree),
    ...         projection="AIT") as wcs:
    ...     ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ...     # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    ...     moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green")
    ...     moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    >>> plt.xlabel('ra')
    >>> plt.ylabel('dec')
    >>> plt.grid(color="black", linestyle="dotted")
    """

    def __init__(
        self,
        fig,
        fov,
        center=coordinates.SkyCoord(0, 0, unit="deg", frame="icrs"),
        coordsys="icrs",
        projection="AIT",
        rotation=coordinates.Angle(0, u.radian),
    ):
        self.w = wcs.WCS(naxis=2)

        width_px, height_px = fig.get_size_inches() * float(fig.dpi)

        cdelt_x = fov.to_value("deg") / float(width_px)

        self.w.wcs.crpix = [width_px / 2.0, height_px / 2.0]
        self.w.wcs.cdelt = [-cdelt_x, cdelt_x]

        if coordsys == "icrs":
            self.w.wcs.crval = [center.icrs.ra.deg, center.icrs.dec.deg]
            self.w.wcs.ctype = ["RA---" + projection, "DEC--" + projection]
        elif coordsys == "galactic":
            self.w.wcs.crval = [center.galactic.l.deg, center.galactic.b.deg]
            self.w.wcs.ctype = ["GLON-" + projection, "GLAT-" + projection]

        theta = rotation.radian
        self.w.wcs.pc = [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ]

    def __enter__(self):
        return self.w

    def __exit__(self, exception_type, exception_value, traceback):
        pass
