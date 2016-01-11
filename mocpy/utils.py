import math

def radec2thetaphi(ra, dec):
    """
    convert equatorial ra, dec in degrees
    to polar theta, phi in radians
    """
    return math.pi/2 - math.radians(dec), math.radians(ra)

def uniq2orderipix(uniq):
    """
    convert a HEALPix pixel coded as a NUNIQ number
    to a (norder, ipix) tuple
    """
    order =  int( (math.log(uniq//4, 2)) // 2 )
    ipix = uniq - 4 * (4**order)
    
    return (order, ipix) 
