import astropy.units as u
import astropy_healpix as ah
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits

from mocpy import MOC, WCS

fits_image_filename = "./../../resources/bayestar.multiorder.fits"

max_order = None
with fits.open(fits_image_filename) as hdul:
    hdul.info()
    hdul[1].columns

    data = hdul[1].data
    max_order = hdul[1].header["MOCORDER"]

uniq = data["UNIQ"]
probdensity = data["PROBDENSITY"]


level, ipix = ah.uniq_to_level_ipix(uniq)
area = ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.steradian)

prob = probdensity * area


cumul_to = np.linspace(0.5, 0.9, 5)[::-1]
colors = ["blue", "green", "yellow", "orange", "red"]
mocs = [MOC.from_valued_healpix_cells(uniq, prob, max_order, cumul_to=c) for c in cumul_to]


# Plot the MOC using matplotlib


fig = plt.figure(111, figsize=(15, 10))
# Define a astropy WCS easily
with WCS(
    fig,
    fov=50 * u.deg,
    center=SkyCoord(315, 15, unit="deg", frame="icrs"),
    coordsys="icrs",
    rotation=Angle(0, u.degree),
    projection="AIT",
) as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    for (moc, c, col) in zip(mocs, cumul_to, colors):
        moc.fill(
            ax=ax,
            wcs=wcs,
            alpha=0.5,
            linewidth=0,
            fill=True,
            color=col,
            label="confidence probability " + str(round(c * 100)) + "%",
        )
        moc.border(ax=ax, wcs=wcs, alpha=0.5, color=col)

    ax.legend()

plt.xlabel("ra")
plt.ylabel("dec")
plt.title("Bayestar")
plt.grid(color="black", linestyle="dotted")
plt.show()
