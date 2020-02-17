from mocpy import MOC, World2ScreenMPL
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
# Load Galex and SDSS
sdss = MOC.from_fits('./../../resources/P-SDSS9-r.fits')
galex = MOC.from_fits('./../../resources/P-GALEXGR6-AIS-FUV.fits')
# Compute their intersection
inter = sdss.intersection(galex)
union = sdss.union(galex)
# Plot the MOC using matplotlib
import matplotlib.pyplot as plt
fig = plt.figure(111, figsize=(10, 10))
# Define a astropy WCS easily
with World2ScreenMPL(fig, 
        fov=160 * u.deg,
        center=SkyCoord(0, 0, unit='deg', frame='icrs'),
        coordsys="icrs",
        rotation=Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    union.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="red", linewidth=0, label="Union")
    union.border(ax=ax, wcs=wcs, alpha=1, color="red")

    inter.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green", linewidth=0, label="Intersection")
    inter.border(ax=ax, wcs=wcs, alpha=1, color="green")
    ax.legend()

plt.xlabel('ra')
plt.ylabel('dec')
plt.title('Logical operations between SDSS and GALEX')
plt.grid(color="black", linestyle="dotted")
plt.show()