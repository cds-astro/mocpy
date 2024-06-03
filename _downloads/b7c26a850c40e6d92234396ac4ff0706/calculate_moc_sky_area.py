import numpy as np
import astropy.units as u
from mocpy import MOC

ALL_SKY = 4 * np.pi * u.steradian
moc = MOC.from_string("0/0-11") # any other MOC instance will work
area = moc.sky_fraction * ALL_SKY
area = area.to(u.deg**2).value
