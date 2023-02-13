import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import Angle, SkyCoord

from mocpy import MOC, WCS

# Load a MOC
filename = "./../../resources/P-SDSS9-r.fits"
moc = MOC.from_fits(filename)
# Plot the MOC using matplotlib
fig = plt.figure(111, figsize=(15, 10))
# Define a astropy WCS easily
wcs = moc.wcs(fig, coordsys="icrs", rotation=Angle(0, u.degree), projection="AIT")
ax = fig.add_subplot(1, 1, 1, projection=wcs)
# Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green")
moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
plt.xlabel("ra")
plt.ylabel("dec")
plt.title("Coverage of P-SDSS9-r")
plt.grid(color="black", linestyle="dotted")
plt.show()
