import pytest
import tempfile
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import ICRS
from astroquery.vizier import Vizier
from astropy.io import fits
from astropy_healpix import HEALPix
import copy

from ..moc import MOC

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

    moc = MOC.from_lonlat(lon=lon, lat=lat, max_norder=7)


def test_moc_from_fits():
    fits_path = 'notebooks/demo-data/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.from_fits(fits_path)



def test_moc_from_fits_images():
    image_path = 'notebooks/demo-data/image_with_mask.fits.gz'

    moc = MOC.from_fits_images([image_path],
                                max_norder=10)


@pytest.fixture()
def moc_from_fits_image():
    image_path = 'notebooks/demo-data/image_with_mask.fits.gz'

    with fits.open(image_path) as hdulist:
        moc = MOC.from_image(header=hdulist[0].header,
                             max_norder=7,
                             mask_arr=hdulist[0].data)
    return moc


def test_moc_from_fits_image(moc_from_fits_image):
    assert moc_from_fits_image


def test_moc_write_and_from_json(moc_from_fits_image):
    tmp_file = tempfile.NamedTemporaryFile()
    moc_from_fits_image.write(tmp_file.name, format='json', write_to_file=True)

    with open(tmp_file.name, 'r') as moc_file:
        import json
        moc_d = json.load(moc_file)
        moc2 = MOC.from_json(json_moc=moc_d)
        assert moc_from_fits_image == moc2


def test_moc_write_to_fits(moc_from_fits_image):
    hdulist = moc_from_fits_image.write(format='fits')
    assert isinstance(hdulist, fits.hdu.hdulist.HDUList)


def test_moc_write_to_json(moc_from_fits_image):
    moc_json = moc_from_fits_image.write(format='json')
    assert isinstance(moc_json, dict)


def test_moc_contains():
    order = 4
    size = 20
    healpix_arr = np.random.randint(0, 12*4**order, size)
    all_healpix_arr = np.arange(12*4**order)
    healpix_outside_arr = np.setdiff1d(all_healpix_arr, healpix_arr)

    moc = MOC.from_json(json_moc={str(order): list(healpix_arr)})

    #viz = Vizier(columns=['*', '_RAJ2000', '_DEJ2000'])
    #viz.ROW_LIMIT = -1
    #table = viz.get_catalogs('I/293/npm2cros')[0]

    hp = HEALPix(nside=(1 << order), order='nested', frame=ICRS())
    lon, lat = hp.healpix_to_lonlat(healpix_arr)
    lon_out, lat_out = hp.healpix_to_lonlat(healpix_outside_arr)

    should_be_inside_arr = moc.contains(ra=lon, dec=lat)
    assert should_be_inside_arr.all()
    should_be_outside_arr = moc.contains(ra=lon_out, dec=lat_out)
    assert not should_be_outside_arr.any()

    # test keep_inside field
    should_be_outside_arr = moc.contains(ra=lon, dec=lat, keep_inside=False)
    assert not should_be_outside_arr.any()
    should_be_inside_arr = moc.contains(ra=lon_out, dec=lat_out, keep_inside=False)
    assert should_be_inside_arr.all()


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


def test_moc_skyfraction():
    moc = MOC.from_json({
        '0': [0, 1, 2, 3, 4, 5]
    })
    assert moc.sky_fraction == 0.5


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


def test_moc_complement():
    moc = MOC.from_fits('notebooks/demo-data/P-GALEXGR6-AIS-FUV.fits')
    assert moc.complement().complement() == moc
