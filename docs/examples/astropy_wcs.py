"""An example for plotting the union of two MOCs with the astropy.wcs module."""

import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
import healpy as hp
from mocpy import MOC
import astropy.units as u

# Generate an abritrary WCS object
wcs = WCS(naxis=2)
wcs.wcs.ctype = ["RA---AIT", "DEC--AIT"]
wcs.wcs.crval = [0.0, 0.0]
wcs.wcs.cdelt = np.array([-0.675, 0.675])
wcs.wcs.crpix = [240.5, 120.5]

# Get HEALPix cell centers
nside = 64
ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest=False, lonlat=True)

# MPL figure with WCS projection
fig = plt.figure(figsize=(15, 7))
ax = plt.subplot(projection=wcs)

# Plot the HEALPix centers (just to test)
ax.scatter(ra, dec, 0.1, transform=ax.get_transform("world"))
ax.grid()

# Define a MOC union
moc1 = MOC.from_cone(
    lon=60 * u.deg, lat=60 * u.deg, radius=20 * u.deg, max_depth=10, delta_depth=2
)
moc2 = MOC.from_cone(
    lon=0 * u.deg, lat=0 * u.deg, radius=10 * u.deg, max_depth=10, delta_depth=2
)
moc = moc1.union(moc2)

# Plot the MOC in the current figure
moc.fill(ax=ax, wcs=wcs)
