import pytest

from ..stmoc import STMOC

import numpy as np

def test_serialization():
    decals = STMOC.from_fits('resources/STMOC/STMoc-DECaLS-g.fits')
    decals_bis = STMOC.load('resources/STMOC/STMoc-DECaLS-g.fits', format='fits')
    assert(decals == decals_bis)

    # Serialize to FITS
    hdulist = decals.serialize(format="fits", pre_v2=True)
    # Deserialize from FITS
    decals_result = decals.deserialization(hdulist)
    assert(decals == decals_result)

    # Save to FITS
    decals.save(path='resources/STMOC/STMoc-DECaLS-g.v2.fits', format='fits', overwrite=True)
    # Load from FITS
    decals_result = STMOC.load(path='resources/STMOC/STMoc-DECaLS-g.v2.fits', format='fits')

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
    assert(decals.max_depth == (61, 9))


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


#### TESTING new features ####
def test_stmoc_save_load_deser():
    stmoc = STMOC.from_string("t61/1 3 5 s3/1-3 t61/50 52 s4/25", 'ascii');
    stmoc_ascii = stmoc.to_string('ascii')
    stmoc_ascii
    stmoc_json = stmoc.to_string('json')
    stmoc_json
    stmoc_bis = STMOC.from_string(stmoc_json, 'json')
    assert stmoc == stmoc_bis
    
    stmoc_bis = STMOC.load('resources/MOC2.0/stmoc.ascii.txt', 'ascii')
    assert stmoc == stmoc_bis
    
    stmoc_bis = STMOC.load('resources/MOC2.0/STMOC.fits', 'fits')
    assert stmoc == stmoc_bis
    
    stmoc.save('resources/MOC2.0/stmoc.py.test.fits', format='fits', overwrite=True)
    stmoc.save('resources/MOC2.0/stmoc.py.test.json', format='json', overwrite=True)
    stmoc.save('resources/MOC2.0/stmoc.py.test.ascii', format='ascii', overwrite=True)
    stmoc_bis = STMOC.load('resources/MOC2.0/stmoc.py.test.fits', 'fits')
    assert stmoc == stmoc_bis
    stmoc_bis = STMOC.load('resources/MOC2.0/stmoc.py.test.json', 'json')
    assert stmoc == stmoc_bis
    stmoc_bis = STMOC.load('resources/MOC2.0/stmoc.py.test.ascii', 'ascii')
    assert stmoc == stmoc_bis
