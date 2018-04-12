import numpy as np
import math

def radec2thetaphi(ra, dec):
    """
    convert equatorial ra, dec in degrees
    to polar theta, phi in radians
    """
    return np.pi/2 - np.radians(dec), np.radians(ra)


def uniq2orderipix(uniq):
    """
    convert a HEALPix pixel coded as a NUNIQ number
    to a (norder, ipix) tuple
    """
    order = int((math.log(uniq//4, 2)) // 2)
    ipix = uniq - 4 * (4**order)
    
    return order, ipix


def orderipix2uniq(n_order, n_pix):
    return n_pix + ((4**n_order) << 2)
