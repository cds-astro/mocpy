import matplotlib.pyplot as plt
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.wcs import WCS
from mocpy import MOC

# Load a MOC
filename = "./../../resources/P-SDSS9-r.fits"
moc = MOC.from_fits(filename)
# Plot the MOC using matplotlib
fig = plt.figure(111, figsize=(15, 10))
# Define a astropy WCS
wcs = WCS(
    {
        "naxis": 2,
        "naxis1": 1620,
        "naxis2": 810,
        "crpix1": 810.5,
        "crpix2": 405.5,
        "cdelt1": -0.2,
        "cdelt2": 0.2,
        "ctype1": "RA---AIT",
        "ctype2": "DEC--AIT",
    },
)
ax = fig.add_subplot(1, 1, 1, projection=wcs, frame_class=EllipticalFrame)
# Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green")
moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
ax.set_aspect(1.0)
plt.xlabel("ra")
plt.ylabel("dec")
plt.title("Coverage of P-SDSS9-r")
plt.grid(color="black", linestyle="dotted")
plt.show()
