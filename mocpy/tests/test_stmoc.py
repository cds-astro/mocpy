import pytest

from ..stmoc import STMOC

import numpy as np

def test_serialization():
    decals = STMOC.from_fits('resources/STMOC/STMoc-DECaLS-g.fits')
    # Serialize to FITS
    hdulist = decals.serialize(format="fits")

    # Deserialize from FITS
    decals_result = decals.deserialization(hdulist)

    assert(decals == decals_result)


from astropy.time import Time
import astropy.units as u
def test_from_times_lonlat():
    times = Time([2440587.50000], format='mjd', scale='tdb')
    lon = [0] * u.deg
    lat = [0] * u.deg

    stmoc = STMOC.from_times_positions(times, 2, lon, lat, 0)

    assert(stmoc.contains(times, lon, lat).all())
    assert(stmoc.contains(times, [180] * u.deg, [0] * u.deg, inside=False).all())
    assert(not stmoc.is_empty())


def test_is_empty():
    empty_stmoc = STMOC()
    assert(empty_stmoc.is_empty())


def test_max_depth():
    decals = STMOC.from_fits('resources/STMOC/STMoc-DECaLS-g.fits')
    assert(decals.max_depth == (29, 9))


def test_union_decals():
    decals = STMOC.from_fits('resources/STMOC/STMoc-DECaLS-g.fits')
    
    result = decals.union(decals)

    assert(decals == result)

def test_intersection_decals():
    decals = STMOC.from_fits('resources/STMOC/STMoc-DECaLS-g.fits')
    
    result = decals.intersection(decals)

    assert(decals == result)

def test_difference_decals():
    decals = STMOC.from_fits('resources/STMOC/STMoc-DECaLS-g.fits')
    
    result = decals.difference(decals)
    assert(result == STMOC())
