from ..stmoc import STMOC
from ..tmoc import TimeMOC
from ..moc import MOC
from astropy.time import Time, TimeDelta
from astropy.table import Table
import astropy.units as u

import pytest
import cdshealpix
import numpy as np

@pytest.fixture()
def stmoc_2mass():
    two_mass_data = Table.read("resources/STMOC/2MASS-list-images.fits.gz", format="fits")

    # Definition of the times, longitudes and latitudes
    times_2mass = Time(two_mass_data["mjd"].data, format="mjd", scale="tdb")
    lon_2mass = two_mass_data["ra"].quantity
    lat_2mass = two_mass_data["dec"].quantity

    time_depth = 10
    spatial_depth = 7
    return STMOC.from_times_positions(
        times_2mass, time_depth, lon_2mass, lat_2mass, spatial_depth
    )


@pytest.fixture()
def stmoc_xmm_dr8():
    xmm_dr8_data = Table.read("resources/STMOC/vizier_votable.b64")
    times_xmm = Time(xmm_dr8_data["MJD0"].data, format="mjd", scale="tdb")
    lon_xmm = xmm_dr8_data["RA_ICRS"].quantity
    lat_xmm = xmm_dr8_data["DE_ICRS"].quantity

    # Create the STMOC
    time_depth = 10
    spatial_depth = 7
    return STMOC.from_times_positions(
        times_xmm, time_depth, lon_xmm, lat_xmm, spatial_depth
    )


def test_serialization():
    decals = STMOC.from_fits("resources/STMOC/STMoc-DECaLS-g.fits")
    decals_bis = STMOC.load("resources/STMOC/STMoc-DECaLS-g.fits", format="fits")
    assert decals == decals_bis

    # Save to FITS
    decals.save(
        path="resources/STMOC/STMoc-DECaLS-g.v2.fits", format="fits", overwrite=True
    )
    # Load from FITS
    decals_result = STMOC.load(
        path="resources/STMOC/STMoc-DECaLS-g.v2.fits", format="fits"
    )

    assert decals == decals_result


def test_from_times_lonlat():
    times = Time([2440587.50000], format="mjd", scale="tdb")
    lon = [0] * u.deg
    lat = [0] * u.deg

    stmoc = STMOC.from_times_positions(times, 2, lon, lat, 0)

    assert stmoc.contains(times, lon, lat).all()
    assert stmoc.contains(times, [180] * u.deg, [0] * u.deg, inside=False).all()
    assert not stmoc.is_empty()


def test_max_depth():
    decals = STMOC.from_fits("resources/STMOC/STMoc-DECaLS-g.fits")
    assert decals.max_depth == (28, 9)


def test_union_decals():
    decals = STMOC.from_fits("resources/STMOC/STMoc-DECaLS-g.fits")

    result = decals.union(decals)

    assert decals == result


def test_intersection_decals():
    decals = STMOC.from_fits("resources/STMOC/STMoc-DECaLS-g.fits")

    result = decals.intersection(decals)

    assert decals == result


def test_difference_decals():
    decals = STMOC.from_fits("resources/STMOC/STMoc-DECaLS-g.fits")

    result = decals.difference(decals)
    assert result.is_empty()


def test_min_max_times(stmoc_xmm_dr8):
    max_time = stmoc_xmm_dr8.max_time
    min_time = stmoc_xmm_dr8.min_time

    tmoc1 = TimeMOC.from_time_ranges(
        min_times = min_time,
        max_times = max_time,
        delta_t = TimeDelta(1e-6, scale="tdb", format="sec")
    )

    one_day = TimeDelta(100, scale="tdb", format="jd");
    tmoc2 = TimeMOC.from_time_ranges(
        min_times = min_time - one_day,
        max_times = max_time + one_day,
        delta_t = TimeDelta(1e-6, scale="tdb", format="sec")
    )

    smoc1 = stmoc_xmm_dr8.query_by_time(tmoc1)
    smoc2 = stmoc_xmm_dr8.query_by_time(tmoc2)
    assert smoc1 == smoc2

def test_query_time(stmoc_xmm_dr8):
    smoc = stmoc_xmm_dr8.query_by_time(TimeMOC.new_empty(max_depth=12))
    assert smoc.empty()

    min_time = stmoc_xmm_dr8.min_time
    tmoc = TimeMOC.from_time_ranges(
        min_times = min_time - TimeDelta(100, scale="tdb", format="jd"),
        max_times = min_time,
        delta_t = TimeDelta(1e-6, scale="tdb", format="sec")
    )

    smoc = stmoc_xmm_dr8.query_by_time(tmoc)
    assert smoc.empty()


def test_stmoc_from_url():
    decals_url = STMOC.from_fits("http://aladin.u-strasbg.fr/java/stmoc/STMoc-DECaLS-g.fits")
    assert decals_url.max_depth == (28, 9)

    decals_local = STMOC.from_fits("resources/STMOC/STMoc-DECaLS-g.fits")
    assert decals_url == decals_local


def test_stmoc_from_time_ranges_positions():
    times_start = Time([2 / TimeMOC.DAY_MICRO_SEC, 3 / TimeMOC.DAY_MICRO_SEC], format="jd", scale="tdb")
    times_end = Time([3 / TimeMOC.DAY_MICRO_SEC, 9 / TimeMOC.DAY_MICRO_SEC], format="jd", scale="tdb")

    time_depth = 61
    spatial_depth = 29

    lon, lat = cdshealpix.healpix_to_lonlat(ipix=np.array([0, 0]), depth=np.array([29, 29]))
    stmoc = STMOC.from_time_ranges_positions(times_start, times_end, lon, lat, time_depth, spatial_depth)

    expected_stmoc = STMOC.from_string("t59/1 60/1 61/8 s29/0")
    assert stmoc == expected_stmoc


def test_stmoc_from_spatial_coverages():
    times_start = Time([2 / TimeMOC.DAY_MICRO_SEC, 3 / TimeMOC.DAY_MICRO_SEC], format="jd", scale="tdb")
    times_end = Time([3 / TimeMOC.DAY_MICRO_SEC, 9 / TimeMOC.DAY_MICRO_SEC], format="jd", scale="tdb")

    time_depth = 61

    smoc = MOC.from_json({"28": [0]})

    stmoc = STMOC.from_spatial_coverages(times_start, times_end, [smoc, smoc], time_depth=time_depth)

    expected_stmoc = STMOC.from_string("t59/1 60/1 61/8 s28/0")
    assert stmoc == expected_stmoc


# --- TESTING new features ---#
def test_stmoc_save_load_deser():
    stmoc = STMOC.from_string("t61/1 3 5 s3/1-3 4/ t61/50 52 s4/25", "ascii")
    stmoc_ascii = stmoc.to_string("ascii")
    stmoc_ascii
    stmoc_json = stmoc.to_string("json")
    stmoc_json
    stmoc_bis = STMOC.from_string(stmoc_json, "json")
    assert stmoc == stmoc_bis

    stmoc_bis = STMOC.load("resources/MOC2.0/stmoc.ascii.txt", "ascii")
    assert stmoc == stmoc_bis

    stmoc_bis = STMOC.load("resources/MOC2.0/STMOC.fits", "fits")
    assert stmoc == stmoc_bis

    stmoc.save("resources/MOC2.0/stmoc.py.test.fits", format="fits", overwrite=True)
    stmoc.save("resources/MOC2.0/stmoc.py.test.json", format="json", overwrite=True)
    stmoc.save("resources/MOC2.0/stmoc.py.test.ascii", format="ascii", overwrite=True)
    stmoc_bis = STMOC.load("resources/MOC2.0/stmoc.py.test.fits", "fits")
    assert stmoc == stmoc_bis
    stmoc_bis = STMOC.load("resources/MOC2.0/stmoc.py.test.json", "json")
    assert stmoc == stmoc_bis
    stmoc_bis = STMOC.load("resources/MOC2.0/stmoc.py.test.ascii", "ascii")
    assert stmoc == stmoc_bis
