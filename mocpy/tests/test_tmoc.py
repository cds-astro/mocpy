import pytest

from astropy.time import Time
from astropy.io import ascii
from ..tmoc import TimeMOC

import numpy as np

def test_empty_tmoc():
    times = Time([], format='jd', scale='tdb')
    tmoc = TimeMOC.from_times(times)
    assert tmoc.empty()
    assert tmoc.total_duration == 0

    with pytest.raises(ValueError):
        min_time = tmoc.min_time

    with pytest.raises(ValueError):
        max_time = tmoc.max_time

    tmoc_ranges = TimeMOC.from_time_ranges(times, times)
    assert tmoc_ranges.empty()
    assert tmoc_ranges.total_duration == 0

def test_simple_tmoc():
    times = Time([2/TimeMOC.DAY_MICRO_SEC, 7/TimeMOC.DAY_MICRO_SEC], format='jd', scale='tdb')
    tmoc = TimeMOC.from_times(times, delta_t=TimeMOC.order_to_time_resolution(29))
    assert tmoc.total_duration.sec == 2 * 1e-6
    assert tmoc.max_order == 29

    tmoc.write('tmoc.txt', format='json', overwrite=True)

def test_tmoc_from_time_ranges():
    """
    Assert a correct tmoc loaded from a fits file is equal to the tmoc built from a CSV file
    containing a list of time intervals
    """
    tmoc = TimeMOC.from_fits('resources/TMOC/HST_SDSSg/TMoc.fits')

    # Load HST_SDSSg from a CSV file
    data = ascii.read('resources/TMOC/HST_SDSSg/uniq-times.csv', format='csv')
    tmoc2 = TimeMOC.from_time_ranges(Time(data['t_min'], format="mjd", scale="tdb"),
                                     Time(data['t_max'], format="mjd", scale="tdb"),
                                     delta_t=TimeMOC.order_to_time_resolution(29))

    assert tmoc == tmoc2

    with pytest.raises(ValueError):
        tmoc_ranges = TimeMOC.from_time_ranges(
            Time([], format="jd", scale="tdb"),
            Time([3], format="jd", scale="tdb")
        )

def test_add_neighbours():
    times = Time([2/TimeMOC.DAY_MICRO_SEC, 7/TimeMOC.DAY_MICRO_SEC], format='jd', scale='tdb')
    times_expected = Time(np.array([1, 2, 3, 6, 7, 8])/TimeMOC.DAY_MICRO_SEC, format='jd', scale='tdb')
    tmoc = TimeMOC.from_times(times, delta_t=TimeMOC.order_to_time_resolution(29))
    
    tmoc_expected = TimeMOC.from_times(times_expected, delta_t=TimeMOC.order_to_time_resolution(29))
    tmoc.add_neighbours()

    assert tmoc == tmoc_expected

def test_remove_neighbours():
    times = Time(np.array([1, 2, 3, 6, 7, 8])/TimeMOC.DAY_MICRO_SEC, format='jd', scale='tdb')
    times_expected = Time([2/TimeMOC.DAY_MICRO_SEC, 7/TimeMOC.DAY_MICRO_SEC], format='jd', scale='tdb')

    tmoc = TimeMOC.from_times(times, delta_t=TimeMOC.order_to_time_resolution(29))
    tmoc_expected = TimeMOC.from_times(times_expected, delta_t=TimeMOC.order_to_time_resolution(29))

    tmoc.remove_neighbours()

    assert tmoc == tmoc_expected

def test_add_remove_back_and_forth():
    times = Time([2/TimeMOC.DAY_MICRO_SEC, 7/TimeMOC.DAY_MICRO_SEC], format='jd', scale='tdb')

    tmoc = TimeMOC.from_times(times, delta_t=TimeMOC.order_to_time_resolution(29))
    tmoc_expected = TimeMOC.from_times(times, delta_t=TimeMOC.order_to_time_resolution(29))

    tmoc.add_neighbours().remove_neighbours()

    assert tmoc == tmoc_expected

def test_contains():
    tmoc = TimeMOC.from_time_ranges(Time(np.array([0]), format="mjd", scale="tdb"),
                                    Time(np.array([1]), format="mjd", scale="tdb"),
                                    delta_t=TimeMOC.order_to_time_resolution(29))

    times_inside = Time(np.linspace(0, 1, num=100), format='mjd', scale='tdb')
    times_outside = Time(np.linspace(1.01, 2, num=100), format='mjd', scale='tdb')
    times_in_and_out = Time(np.linspace(0.9, 2, num=100), format='mjd', scale='tdb')

    assert tmoc.contains(times_inside).all()
    assert (~tmoc.contains(times_outside)).all()
    assert tmoc.contains(times_in_and_out).any()

@pytest.mark.parametrize("a, b, expect", [
    (TimeMOC.from_times(Time(np.array([0, 2]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([0, 1, 2]), format='jd', scale='tdb'))),

    (TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb'))),

    (TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb'))),
])
def test_union(a, b, expect):
    res = a.union(b)

    assert res == expect

@pytest.mark.parametrize("a, b, expect", [
    (TimeMOC.from_times(Time(np.array([0, 2]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([0, 2]), format='jd', scale='tdb'))),

    (TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb'))),

    (TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb'))),
])
def test_difference(a, b, expect):
    res = a.difference(b)

    assert res == expect

@pytest.mark.parametrize("a, b, expect", [
    (TimeMOC.from_times(Time(np.array([0, 1, 2]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb'))),

    (TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([1]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb'))),

    (TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb')),
     TimeMOC.from_times(Time(np.array([]), format='jd', scale='tdb'))),
])
def test_intersection(a, b, expect):
    res = a.intersection(b)

    assert res == expect
