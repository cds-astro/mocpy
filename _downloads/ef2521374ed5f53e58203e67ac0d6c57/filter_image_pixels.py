import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.wcs import WCS as WCS

from mocpy import MOC

# load 2MASS cutout covering the galactic plane
hdu = fits.open(
    "http://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips=CDS%2FP%2F2MASS%2FK&width=1200&height=700&fov=30&projection=TAN&coordsys=galactic&rotation_angle=0.0&object=gal%20center&format=fits"
)

# load Spitzer MOC
moc = MOC.from_fits("http://skies.esac.esa.int/Spitzer/IRAC1_bright_ISM/Moc.fits")

# create WCS from 2MASS image header
twomass_wcs = WCS(header=hdu[0].header)

# compute skycoords for every pixel of the image
width = hdu[0].header["NAXIS1"]
height = hdu[0].header["NAXIS2"]

xv, yv = np.meshgrid(np.arange(0, width), np.arange(0, height))
skycoords = twomass_wcs.pixel_to_world(xv, yv)

ra, dec = skycoords.icrs.ra.deg, skycoords.icrs.dec.deg
# Get the skycoord lying in the MOC mask
mask_in_moc = moc.contains(ra * u.deg, dec * u.deg)

img = hdu[0].data
img_test = img.copy()

# Set to zero the pixels that do not belong to the MOC
img_test[~mask_in_moc] = 0

# Plot the filtered pixels image
fig = plt.figure(111, figsize=(15, 10))
ax = fig.add_subplot(1, 1, 1, projection=twomass_wcs)

im = ax.imshow(
    img_test, origin="lower", norm=simple_norm(hdu[0].data, "sqrt", min_cut=-1, max_cut=150),
)
plt.show()

# Sum of pixels lying in the MOC
np.sum(img[mask_in_moc])

# Sum of pixels of the whole image
np.sum(img)
