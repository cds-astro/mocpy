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
    order = ((np.log2(uniq//4)) // 2)
    order = order.astype(int)
    ipix = uniq - 4 * (4**order)
    
    return order, ipix


def orderipix2uniq(n_order, n_pix):
    return n_pix + ((4**n_order) << 2)


def number_trailing_zeros(i):
    # O(log(i)) complexity here. Can be done in O(1) using more low-level algo
    nb = 0
    while (i & 1) == 0 and i > 0:
        nb += 1
        i = i >> 1

    return nb

