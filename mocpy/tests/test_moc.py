import pytest
import copy
import sys

import numpy as np

from astropy.coordinates import SkyCoord, ICRS, Angle
from astropy.io.votable import parse_single_table
import astropy.units as u
from astropy.io import fits

import cdshealpix

from ..moc import MOC, WCS


#### TESTING MOC creation ####
def get_random_skycoords(size):
    return SkyCoord(ra=np.random.uniform(0, 360, size),
                    dec=np.random.uniform(-90, 90, size),
                    unit="deg")


skycoords1 = get_random_skycoords(size=1000)
skycoords2 = get_random_skycoords(size=2000)
skycoords3 = get_random_skycoords(size=50000)

@pytest.fixture()
def skycoords_gen_f():
    def gen_f(size):
        return SkyCoord(np.random.uniform(0, 360, size), np.random.uniform(-90, 90, size), unit='deg')

    return gen_f


@pytest.fixture()
def lonlat_gen_f():
    def gen_f(size):
        return np.random.uniform(0, 360, size) * u.deg, np.random.uniform(-90, 90, size) * u.deg

    return gen_f


@pytest.mark.parametrize("size", [
    1000,
    10000,
    50000
])
def test_moc_from_skycoords(skycoords_gen_f, size):
    skycoords = skycoords_gen_f(size)
    moc = MOC.from_skycoords(skycoords, max_norder=7)


@pytest.mark.parametrize("size", [
    1000,
    10000,
    50000
])
def test_moc_from_lonlat(lonlat_gen_f, size):
    lon, lat = lonlat_gen_f(size)
    moc = MOC.from_lonlat(lon=lon, lat=lat, max_norder=6)


def test_from_healpix_cells():
    ipix = np.array([40, 87, 65])
    depth = np.array([3, 3, 3])
    fully_covered = np.array([True, True, True])

    moc = MOC.from_healpix_cells(ipix, depth, fully_covered)


def test_moc_from_fits():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.from_fits(fits_path)


def test_moc_consistent_with_aladin():
    truth = MOC.from_fits('resources/CDS-I-125A-catalog_MOC.fits')
    table = parse_single_table("resources/I_125A_catalog.vot").to_table()

    moc = MOC.from_lonlat(
        table['_RAJ2000'].T * u.deg,
        table['_DEJ2000'].T * u.deg,
        max_norder=8
    )

    assert moc == truth


def test_moc_from_fits_images():
    image_path = 'resources/image_with_mask.fits.gz'

    moc = MOC.from_fits_images([image_path],
                                max_norder=10)


@pytest.fixture()
def moc_from_fits_image():
    image_path = 'resources/image_with_mask.fits.gz'

    with fits.open(image_path) as hdulist:
        moc = MOC.from_image(header=hdulist[0].header,
                             max_norder=7,
                             mask=hdulist[0].data)
    return moc


@pytest.fixture()
def moc_from_json():
    return MOC.from_json({'8': [45, 78], '4': [42, 57]})


def test_moc_from_fits_image(moc_from_fits_image):
    assert isinstance(moc_from_fits_image, MOC)


def test_moc_serialize_and_from_json(moc_from_json):
    ipix_d = moc_from_json.serialize(format="json")
    moc2 = MOC.from_json(ipix_d)
    assert moc_from_json == moc2


@pytest.mark.parametrize("expected, moc_str", [
    (MOC.from_json({'5': [8, 9, 10, 42, 43, 44, 45, 54, 46], '6':[4500], '7':[], '8':[45]}),
    '5/8-10,42-46,54,8 6/4500 8/45'),
    (MOC.from_json({}), '0/'),
    (MOC.from_json({'29': [101]}), '29/101'),
    (MOC.from_json({'0': [1, 0, 9]}), '0/0-1,9'),
    (MOC.from_json({'0': [2, 9], '1': [9]}), '0/2,9'),
    (MOC.from_json({'0': [2], '8': [8, 9, 10], '11': []}), '0/2\r ,\n 8/8-10\n 11/'),
])
def test_from_str(expected, moc_str):
    assert MOC.from_str(moc_str) == expected


def test_moc_skyfraction():
    moc = MOC.from_json({
        '0': [0, 1, 2, 3, 4, 5]
    })
    assert moc.sky_fraction == 0.5


#### TESTING MOC serialization ####
def test_moc_serialize_to_fits(moc_from_fits_image):
    hdulist = moc_from_fits_image.serialize(format='fits')
    assert isinstance(hdulist, fits.hdu.hdulist.HDUList)


def test_moc_serialize_to_json(moc_from_fits_image):
    moc_json = moc_from_fits_image.serialize(format='json')
    assert isinstance(moc_json, dict)


@pytest.mark.parametrize("moc, expected", [
    (MOC.from_json({'5': [8, 9, 10, 42, 43, 44, 45, 54, 46], '6':[4500], '7':[], '8':[45]}),
    '5/8-10,42-46,54 6/4500 8/45'),
    (MOC.from_json({}), ''),
    (MOC.from_json({'29': [101]}), '29/101'),
    (MOC.from_json({'0': [1, 0, 9]}), '0/0-1,9'),
    (MOC.from_json({'0': [2, 9], '1': [9]}), '0/2,9'),
])
def test_serialize_to_str(moc, expected):
    assert moc.serialize(format="str") == expected


@pytest.mark.parametrize("filename, overwrite, format, os_error", [
    ('moc', True, 'fits', False),
    ('moc', False, 'fits', True),
    ('moc', True, 'json', False),
    ('moc', True, 'str', False),
    ('moc', False, 'str', True),
])
def test_write(moc_from_json, filename, overwrite, format, os_error):
    if os_error:
        with pytest.raises(OSError):
            moc_from_json.write(filename, format=format, overwrite=overwrite)
    else:
        moc_from_json.write(filename, format=format, overwrite=overwrite)


#### TESTING MOC plot functions ####
def test_mpl_fill():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.from_fits(fits_path)

    import matplotlib.pyplot as plt
    fig = plt.figure(111, figsize=(10, 10))
    with WCS(fig,
         fov=50 * u.deg,
         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
         coordsys="icrs",
         rotation=Angle(0, u.degree),
         projection="AIT") as wcs:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc.fill(ax=ax, wcs=wcs, alpha=0.5, color='r')


def test_mpl_border():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.from_fits(fits_path)

    import matplotlib.pyplot as plt
    fig = plt.figure(111, figsize=(10, 10))
    with WCS(fig,
         fov=50 * u.deg,
         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
         coordsys="icrs",
         rotation=Angle(0, u.degree),
         projection="AIT") as wcs:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc.border(ax=ax, wcs=wcs, color='g')


#### TESTING MOC features ####
def test_moc_contains():
    order = 4
    size = 20
    healpix_arr = np.random.randint(0, 12*4**order, size)
    all_healpix_arr = np.arange(12*4**order)
    healpix_outside_arr = np.setdiff1d(all_healpix_arr, healpix_arr)

    moc = MOC.from_json(json_moc={str(order): healpix_arr.tolist()})

    lon, lat = cdshealpix.healpix_to_lonlat(healpix_arr, order)
    lon_out, lat_out = cdshealpix.healpix_to_lonlat(healpix_outside_arr, order)

    should_be_inside_arr = moc.contains(ra=lon, dec=lat)
    assert should_be_inside_arr.all()
    should_be_outside_arr = moc.contains(ra=lon_out, dec=lat_out)
    assert not should_be_outside_arr.any()

    # test keep_inside field
    should_be_outside_arr = moc.contains(ra=lon, dec=lat, keep_inside=False)
    assert not should_be_outside_arr.any()
    should_be_inside_arr = moc.contains(ra=lon_out, dec=lat_out, keep_inside=False)
    assert should_be_inside_arr.all()


# TODO: IMPROVE THE ALGO
'''
def test_boundaries():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.from_fits(fits_path)
    moc = moc.degrade_to_order(6)
    boundaries_l = moc.get_boundaries()
'''


@pytest.fixture()
def mocs():
    moc1 = {'1': [0]}
    moc1_increased = {'0': [0], '1': [17, 19, 22, 23, 35]}
    moc2 = {'1': [30]}
    moc2_increased = {'0': [7], '1': [8, 9, 25, 43, 41]}

    return dict(moc1=MOC.from_json(moc1),
                moc1_increased=MOC.from_json(moc1_increased),
                moc2=MOC.from_json(moc2),
                moc2_increased=MOC.from_json(moc2_increased))


def test_add_neighbours(mocs):
    mocs['moc1'].add_neighbours()
    assert mocs['moc1'] == mocs['moc1_increased']

    mocs['moc2'].add_neighbours()
    assert mocs['moc2'] == mocs['moc2_increased']


def test_remove_neighbours(mocs):
    mocs['moc1_increased'].remove_neighbours()
    mocs['moc2_increased'].remove_neighbours()
    assert mocs['moc1_increased'] == mocs['moc1']
    assert mocs['moc2_increased'] == mocs['moc2']


def test_neighbours(mocs):
    moc1 = copy.deepcopy(mocs['moc1'])
    moc2 = copy.deepcopy(mocs['moc2'])
    moc1.add_neighbours().remove_neighbours()
    moc2.add_neighbours().remove_neighbours()
    assert moc1 == mocs['moc1']
    assert moc2 == mocs['moc2']


#### TESTING MOC operations ####
@pytest.fixture()
def mocs_op():
    moc1 = MOC.from_json({
        '0': [0, 2, 3, 4, 5]
    })
    moc2 = MOC.from_json({
        '0': [0, 1, 7, 4, 3]
    })
    return dict(first=moc1, second=moc2)


def test_moc_union(mocs_op):
    assert mocs_op['first'].union(mocs_op['second']) == MOC.from_json({
        '0': [0, 1, 2, 3, 4, 5, 7]
    })


def test_moc_intersection(mocs_op):
    assert mocs_op['first'].intersection(mocs_op['second']) == MOC.from_json({
        '0': [0, 3, 4]
    })


def test_moc_difference(mocs_op):
    assert mocs_op['first'].difference(mocs_op['second']) == MOC.from_json({
        '0': [2, 5]
    })


def test_moc_complement_consistency():
    moc = MOC.from_fits('resources/P-GALEXGR6-AIS-FUV.fits')
    assert moc.complement().complement() == moc


@pytest.mark.parametrize("input, expected", [
    (MOC.from_json({'0': [1, 3]}), MOC.from_json({'0': [0, 2, 4, 5, 6, 7, 8, 9, 10, 11]}))
])
def test_moc_complement(input, expected):
    assert input.complement() == expected
