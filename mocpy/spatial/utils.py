from astropy import wcs

def make_wcs(crpix, crval, cdelt, ctype):
    world_coord_sys = wcs.WCS(naxis=2)

    world_coord_sys.wcs.crpix = crpix
    world_coord_sys.wcs.cdelt = cdelt
    world_coord_sys.wcs.crval = crval
    world_coord_sys.wcs.ctype = ctype

    return world_coord_sys
