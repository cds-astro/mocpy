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
