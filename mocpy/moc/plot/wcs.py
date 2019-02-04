import numpy as np

from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import wcs

import astropy.units as u

from matplotlib.pyplot import figure

class WCS:
    def __init__(self,
                 fig,
                 fov,
                 center=SkyCoord(0, 0, unit="deg", frame="icrs"),
                 coordsys="icrs",
                 projection="AIT",
                 rotation=Angle(0, u.radian)):
        """
        Create a WCS for vizualizing a MOC in a matplotlib axis.

        Parameters
        ----------
        fig: `~matplotlib.pyplot.figure`
            The matplotlib figure used for plotting the MOC.
        fov: `~astropy.units.Quantity`
            Field of view
        center: `~astropy.coordinates.SkyCoord`, optional
            World coordinates matching with the center of the plot. Default to (0 deg, 0 deg) (in ICRS frame)
        coordsys: str, optional
            Coordinate system. Default to "icrs". Must be in ["icrs", "galactic"]
        projection: str, optional
            World system -> Image system projection type. See http://docs.astropy.org/en/stable/wcs/#supported-projections for
            the projections currently supported in astropy. Default to Aitoff.
        rotation: `~astropy.coordinates.Angle`, optional
            The angle of rotation. Default to no rotation.

        Returns
        -------
        wcs : `~astropy.wcs.WCS`
            The WCS that can be passed to mocpy.MOC.fill/border.
        """
        self.w = wcs.WCS(naxis=2)
        
        width_px, height_px = fig.get_size_inches() * fig.dpi

        cdelt = fov.to_value("deg")/width_px

        self.w.wcs.crpix = [width_px/2, height_px/2]
        self.w.wcs.cdelt = [-cdelt, cdelt]

        if coordsys == 'icrs':
            self.w.wcs.crval = [center.icrs.ra.deg, center.icrs.dec.deg]
            self.w.wcs.ctype = ['RA---' + projection, 'DEC--' + projection]
        elif coordsys == 'galactic':
            self.w.wcs.crval = [center.galactic.l.deg, center.galactic.b.deg]
            self.w.wcs.ctype = ['GLON-' + projection, 'GLAT-' + projection]

        theta = rotation.radian
        self.w.wcs.pc = [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ]

    def __enter__(self):
        return self.w

    def __exit__(self, exception_type, exception_value, traceback):
        pass

