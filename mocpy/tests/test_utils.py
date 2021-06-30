import pytest
import numpy as np
from astropy.time import Time
from .. import utils

def test_time_to_microsec_1():
    # t = Time([['1998-01-01', '1999-01-01']], format="iso", scale="tdb")

    t = Time('1999-01-01T00:00:00.123456789', format='isot', scale='tdb')
    #t = Time('0000-01-01T00:00:00.123456789', format='isot', scale='tdb')
    #print(t.jd * utils.DAY_MICRO_SEC)
    us = utils.times_to_microseconds(t)
    jd = utils.microseconds_to_times(us)
    assert us==211781908800123456
    assert jd.jd==t.jd

def test_time_to_microsec_2():
    t = Time([['1998-01-01', '1999-01-01']], format="iso", scale="tdb")
    us1 = utils.times_to_microseconds(t)
    us2 = np.asarray(t.jd * 86400000000, dtype=np.uint64)
    jd1 = utils.microseconds_to_times(us1)
    jd2 = utils.microseconds_to_times(us2)
    print(us1)
    print(us2)
    assert (us1==us2).all()
    assert (jd1.jd==jd2.jd).all()
