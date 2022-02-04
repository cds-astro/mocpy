import pytest
import copy
import sys

import numpy as np

from astropy.coordinates import SkyCoord, ICRS, Angle
from astropy.io.votable import parse_single_table
import astropy.units as u
from astropy.io import fits

import cdshealpix

from ..interval_set import IntervalSet
from ..moc import MOC, World2ScreenMPL



def test_interval_min_depth():
    big_cells = np.array([[0, 4**29]], dtype=np.uint64)
    itv_result = MOC(IntervalSet(big_cells, make_consistent=False), min_depth=1)

    small_cells = np.array([[0, 4**28], [4**28, 2*4**28], [2*4**28, 3*4**28], [3*4**28, 4**29]], dtype=np.uint64)
    itv_small_cells = MOC(IntervalSet(small_cells, make_consistent=False), make_consistent=False)
    assert itv_result == itv_small_cells

def test_interval_set_complement():
    assert MOC().complement() == MOC(IntervalSet(np.array([[0, 12*4**29]], dtype=np.uint64)))
    assert MOC().complement().complement() == MOC()
    assert MOC(IntervalSet(np.array([[1, 2], [6, 8], [5, 6]], dtype=np.uint64))).complement() == \
        MOC(IntervalSet(np.array([[0, 1], [2, 5], [8, 12*4**29]], dtype=np.uint64)))

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
    moc = MOC.load(fits_path, 'fits')


def test_moc_consistent_with_aladin():
    truth = MOC.load('resources/CDS-I-125A-catalog_MOC.fits', 'fits')
    table = parse_single_table("resources/I_125A_catalog.vot").to_table()

    moc = MOC.from_lonlat(
        table['_RAJ2000'].T * u.deg,
        table['_DEJ2000'].T * u.deg,
        max_norder=8
    )

    assert moc == truth


def test_moc_from_fits_images():
    image_path = 'resources/image_with_mask.fits.gz'

    moc = MOC.from_fits_images([image_path], max_norder=15)


def test_from_fits_images_2():
    MOC.from_fits_images(['resources/u_gal.fits'], max_norder=10)


@pytest.fixture()
def moc_from_fits_image():
    image_path = 'resources/image_with_mask.fits.gz'

    with fits.open(image_path) as hdulist:
        moc = MOC.from_fits_image(hdu=hdulist[0], max_norder=7, mask=hdulist[0].data)

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
    '5/8-10 42-46 54\n\r 8 6/4500 8/45'),
    (MOC.from_json({}), '0/'),
    (MOC.from_json({'29': [101]}), '29/101'),
    (MOC.from_json({'0': [1, 0, 9]}), '0/0-1 9'),
    (MOC.from_json({'0': [2, 9], '1': [9]}), '0/2 9'),
    (MOC.from_json({'0': [2], '8': [8, 9, 10], '11': []}), '0/2\r \n 8/8-10\n 11/'),
])
def test_from_str(expected, moc_str):
    assert MOC.from_str(moc_str) == expected

@pytest.mark.parametrize("expected, moc_str", [
    (MOC.from_json({'5': [8, 9, 10, 42, 43, 44, 45, 54, 46], '6':[4500], '7':[], '8':[45]}),
    '5/8-10 42-46 54\n\r 6/4500 8/45'),
    (MOC.from_json({}), '0/'),
    (MOC.from_json({'29': [101]}), '29/101'),
    (MOC.from_json({'0': [1, 0, 9]}), '0/0-1 9'),
    (MOC.from_json({'0': [2, 9], '1': [9]}), '0/2 9'),
    (MOC.from_json({'0': [2], '8': [8, 9, 10], '11': []}), '0/2\r \n 8/8-10\n 11/'),
])
def test_from_string(expected, moc_str):
    assert MOC.from_string(moc_str, 'ascii') == expected


def test_moc_full_skyfraction():
    moc = MOC.from_json({
        '0': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    })
    assert moc.sky_fraction == 1.0


def test_moc_skyfraction():
    moc = MOC.from_json({
        '0': [0, 1, 2, 3, 4, 5]
    })
    assert moc.sky_fraction == 0.5


def test_sky_fraction_on_empty_coverage():
    moc = MOC()
    assert moc.sky_fraction == 0


#### TESTING MOC serialization ####
def test_moc_serialize_to_fits(moc_from_fits_image):
    hdulist = moc_from_fits_image.serialize(format='fits')
    assert isinstance(hdulist, fits.hdu.hdulist.HDUList)


def test_moc_serialize_to_json(moc_from_fits_image):
    moc_json = moc_from_fits_image.serialize(format='json')
    assert isinstance(moc_json, dict)


@pytest.mark.parametrize("moc, expected", [
    (MOC.from_json({'5': [8, 9, 10, 42, 43, 44, 45, 54, 46], '6':[4500], '7':[], '8':[45]}),
    '5/8-10 42-46 54 6/4500 8/45'),
    (MOC.from_json({}), ''),
    (MOC.from_json({'29': [101]}), '29/101'),
    (MOC.from_json({'0': [1, 0, 9]}), '0/0-1 9'),
    (MOC.from_json({'0': [2, 9], '1': [9]}), '0/2 9'),
])
def test_serialize_to_str(moc, expected):
    assert moc.serialize(format="str") == expected


@pytest.mark.parametrize("filename, overwrite, format, os_error", [
    ('moc', True, 'fits', False),
    ('moc', False, 'fits', True),
    ('moc', True, 'json', False),
    ('moc', True, 'ascii', False),
    ('moc', False, 'ascii', True),
])
def test_write(moc_from_json, filename, overwrite, format, os_error):
    if os_error:
        with pytest.raises(OSError):
            moc_from_json.save(filename, format=format, overwrite=overwrite)
    else:
        moc_from_json.save(filename, format=format, overwrite=overwrite)


#### TESTING MOC plot functions ####
def test_mpl_fill():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.load(fits_path, 'fits')

    import matplotlib.pyplot as plt
    fig = plt.figure(111, figsize=(10, 10))
    with World2ScreenMPL(fig,
         fov=50 * u.deg,
         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
         coordsys="icrs",
         rotation=Angle(0, u.degree),
         projection="AIT") as wcs:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc.fill(ax=ax, wcs=wcs, alpha=0.5, color='r')


def test_mpl_border():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.load(fits_path, 'fits')

    import matplotlib.pyplot as plt
    fig = plt.figure(111, figsize=(10, 10))
    with World2ScreenMPL(fig,
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


def test_degrade_to_order():
    hst_fits_path = 'resources/hst.fits'
    hst_moc = MOC.load(hst_fits_path, 'fits')

    max_depth = hst_moc.max_order

    for order in reversed(range(0, max_depth)):
        hst_moc = hst_moc.degrade_to_order(order)
        assert(hst_moc.sky_fraction <= 1.0)


def test_from_ring():
     moc = MOC.from_ring(
         lon=0 * u.deg,
         lat=0 * u.deg,
         internal_radius=Angle(5, u.deg),
         external_radius=Angle(10, u.deg),
         max_depth=10)

# TODO: IMPROVE THE ALGO
'''
def test_boundaries():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.load(fits_path, 'fits')
    moc = moc.degrade_to_order(6)
    boundaries_l = moc.get_boundaries()
'''
def test_from_elliptical_cone():
    moc = MOC.from_elliptical_cone(
        lon=0 * u.deg,
        lat=0 * u.deg,
        a=Angle(10, u.deg),
        b=Angle(5, u.deg),
        pa=Angle(0, u.deg),
        max_depth=10)


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
    print(mocs['moc1'])
    mocs['moc1'].add_neighbours()
    print(mocs['moc1'])
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
    moc = MOC.load('resources/P-GALEXGR6-AIS-FUV.fits', 'fits')
    assert moc.complement().complement() == moc


def test_from_fits_old():
    moc = MOC.from_fits('resources/V_147_sdss12.moc.fits')
    assert moc.complement().complement() == moc


@pytest.mark.parametrize("input, expected", [
    (MOC.from_json({'0': [1, 3]}), MOC.from_json({'0': [0, 2, 4, 5, 6, 7, 8, 9, 10, 11]}))
])
def test_moc_complement(input, expected):
    assert input.complement() == expected


def test_spatial_res_to_order():
    order = np.arange(14)

    res = MOC.order_to_spatial_resolution(order)
    output = MOC.spatial_resolution_to_order(res)

    assert (order == output).all()


def test_from_valued_healpix_cells_empty():
    uniq = np.array([])
    values = np.array([])

    MOC.from_valued_healpix_cells(uniq, values)


def test_from_valued_healpix_cells_different_sizes():
    uniq = np.array([500])
    values = np.array([])

    with pytest.raises(ValueError):
        MOC.from_valued_healpix_cells(uniq, values)


def test_from_valued_healpix_cells_cumul_from_sup_cumul_to():
    uniq = np.array([500])
    values = np.array([1.0])

    with pytest.raises(ValueError):
        MOC.from_valued_healpix_cells(uniq, values, cumul_from=0.8, cumul_to=-5.0)


@pytest.mark.parametrize("cumul_from, cumul_to", [
    (-5.0, 1.0),
    (np.nan, np.inf),
    (np.nan, np.nan),
    (np.inf, np.nan),
    (-10.0, -5.0)
])
def test_from_valued_healpix_cells_weird_values(cumul_from, cumul_to):
    uniq = np.array([500])
    values = np.array([-1.0])

    MOC.from_valued_healpix_cells(uniq, values, cumul_from=cumul_from, cumul_to=cumul_to)


def test_from_valued_healpix_cells_bayestar():
    from astropy.io import fits
    fits_image_filename = './resources/bayestar.multiorder.fits'

    with fits.open(fits_image_filename) as hdul:
        hdul.info()
        hdul[1].columns

        data = hdul[1].data

    uniq = data['UNIQ']
    probdensity = data['PROBDENSITY']

    import astropy_healpix as ah
    import astropy.units as u

    level, ipix = ah.uniq_to_level_ipix(uniq)
    area = ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.steradian)

    prob = probdensity * area

    cumul_to = np.linspace(0.01, 2.0, num=10)

    for b in cumul_to:
        MOC.from_valued_healpix_cells(uniq, prob, 12, cumul_from=0.0, cumul_to=b)

def test_from_valued_healpix_cells_bayestar_and_split():
   fits_mom_filename = './resources/bayestar.multiorder.fits'
   moc = MOC.from_multiordermap_fits_file(fits_mom_filename, cumul_from=0.0, cumul_to=0.9)
   count = moc.split_count()
   assert count == 2
   mocs = list(moc.split())
   assert len(mocs)== 2
   for moc in mocs:
     assert moc.max_order == 8


#### TESTING new features ####
def test_moc_save_load_deser():
    smoc = MOC.from_string("3/3 10 4/16-18 22 5/19-20 17/222 28/123456789 29/", 'ascii')
    smoc_ascii = smoc.to_string('ascii')
    smoc_ascii
    smoc_json = smoc.to_string('json')
    smoc_json
    smoc_bis = MOC.from_string(smoc_json, 'json')
    assert smoc == smoc_bis

    smoc_bis = MOC.load('resources/MOC2.0/smoc.ascii.txt', 'ascii')
    assert smoc == smoc_bis

    smoc_bis = MOC.load('resources/MOC2.0/SMOC.fits', 'fits')
    assert smoc == smoc_bis

    smoc.save('resources/MOC2.0/smoc.py.test.fits', format='fits', overwrite=True)
    smoc.save('resources/MOC2.0/smoc.py.test.json', format='json', overwrite=True)
    smoc.save('resources/MOC2.0/smoc.py.test.ascii', format='ascii', overwrite=True)
    smoc_bis = MOC.load('resources/MOC2.0/smoc.py.test.fits', 'fits')
    assert smoc == smoc_bis
    smoc_bis = MOC.load('resources/MOC2.0/smoc.py.test.json', 'json')
    assert smoc == smoc_bis
    smoc_bis = MOC.load('resources/MOC2.0/smoc.py.test.ascii', 'ascii')
    assert smoc == smoc_bis


