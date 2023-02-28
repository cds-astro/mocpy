import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from mocpy import MOC, WCS

moc = MOC.from_fits("polygon_moc.fits")
skycoords = moc.get_boundaries()

# Plot the MOC using matplotlib
fig = plt.figure(111, figsize=(10, 10))
# Define a astropy WCS easily
with WCS(
    fig,
    fov=20 * u.deg,
    center=SkyCoord(10, 5, unit="deg", frame="icrs"),
    coordsys="icrs",
    rotation=Angle(0, u.degree),
    # The gnomonic projection transforms great circles into straight lines.
    projection="TAN",
) as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="red", linewidth=1)
    moc.border(ax=ax, wcs=wcs, alpha=1, color="red")

    # Plot the border
    from astropy.wcs.utils import skycoord_to_pixel

    x, y = skycoord_to_pixel(skycoords[0], wcs)
    p = Path(np.vstack((x, y)).T)
    patch = PathPatch(p, color="black", fill=False, alpha=0.75, lw=2)
    ax.add_patch(patch)

plt.xlabel("ra")
plt.ylabel("dec")
plt.title("Get the border(s) of a MOC")
plt.grid(color="black", linestyle="dotted")
plt.show()
