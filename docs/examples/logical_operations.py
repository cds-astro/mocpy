import matplotlib.pyplot as plt
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.wcs import WCS
from mocpy import MOC

# Load Galex and SDSS
sdss = MOC.from_fits("./../../resources/P-SDSS9-r.fits")
galex = MOC.from_fits("./../../resources/P-GALEXGR6-AIS-FUV.fits")
# Compute their intersection
inter = sdss & galex
union = sdss + galex
# Plot the MOC using matplotlib
fig = plt.figure(111, figsize=(10, 8))
# Define a astropy WCS
wcs = WCS(
    {
        "naxis": 2,
        "naxis1": 3240,
        "naxis2": 1620,
        "crpix1": 1620.5,
        "crpix2": 810.5,
        "cdelt1": -0.1,
        "cdelt2": 0.1,
        "ctype1": "RA---AIT",
        "ctype2": "DEC--AIT",
    },
)
ax = fig.add_subplot(1, 1, 1, projection=wcs, frame_class=EllipticalFrame)
# Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
union.fill(
    ax=ax,
    wcs=wcs,
    alpha=0.5,
    fill=True,
    color="blue",
    linewidth=0,
    label="Union",
)
union.border(ax=ax, wcs=wcs, alpha=1, color="black")

inter.fill(
    ax=ax,
    wcs=wcs,
    alpha=0.5,
    fill=True,
    color="green",
    linewidth=0,
    label="Intersection",
)
inter.border(ax=ax, wcs=wcs, alpha=1, color="black")
ax.legend()
ax.set_aspect(1.0)
plt.xlabel("ra")
plt.ylabel("dec")
plt.title("Logical operations between SDSS and GALEX")
plt.grid(color="black", linestyle="dotted")
plt.show()
