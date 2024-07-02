import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.wcs import WCS
from mocpy import MOC

# this can be found in mocpy's repository
# https://github.com/cds-astro/mocpy/blob/master/resources/bayestar.multiorder.fits
fits_image_filename = "./../../resources/bayestar.multiorder.fits"

with fits.open(fits_image_filename) as hdul:
    hdul.info()
    data = hdul[1].data
    max_order = hdul[1].header["MOCORDER"]

uniq = data["UNIQ"]
probdensity = data["PROBDENSITY"]

# let's convert the probability density into a probability
orders = (np.log2(uniq // 4)) // 2
area = 4 * np.pi / np.array([MOC.n_cells(int(order)) for order in orders]) * u.sr
prob = probdensity * area

# now we create the mocs corresponding to different probability thresholds
cumul_to = np.linspace(0.5, 0.9, 5)[::-1]
colors = ["blue", "green", "yellow", "orange", "red"]
mocs = [
    MOC.from_valued_healpix_cells(uniq, prob, max_order, cumul_to=c) for c in cumul_to
]


# Plot the MOC using matplotlib
fig = plt.figure(111)
# Define a astropy WCS
wcs = WCS(
    {
        "naxis": 2,
        "naxis1": 3240,
        "naxis2": 3240,
        "crpix1": 1620.5,
        "crpix2": 1620.5,
        "cdelt1": -0.0353,
        "cdelt2": 0.0353,
        "ctype1": "RA---SIN",
        "ctype2": "DEC--SIN",
    },
)

ax = fig.add_subplot(1, 1, 1, projection=wcs, frame_class=EllipticalFrame)
# Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
for moc, c, col in zip(mocs, cumul_to, colors):
    moc.fill(
        ax=ax,
        wcs=wcs,
        alpha=0.5,
        linewidth=0,
        fill=True,
        color=col,
        label="confidence proba " + str(round(c * 100)) + "%",
    )
    moc.border(ax=ax, wcs=wcs, alpha=0.5, color=col)
ax.legend(loc="upper center", bbox_to_anchor=(0, 1))
ax.set_aspect(1.0)
plt.xlabel("ra")
plt.ylabel("dec")
plt.title("Bayestar")
plt.grid(color="black", linestyle="dotted")
plt.tight_layout()
plt.show()
