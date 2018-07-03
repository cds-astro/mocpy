import pytest
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from .moc import MOC
from .tmoc import TimeMoc


def get_random_skycoords(size):
    return SkyCoord(ra=np.random.uniform(0, 360, size),
                    dec=np.random.uniform(-90, 90, size),
                    unit="deg")


skycoords1 = get_random_skycoords(size=1000)
skycoords2 = get_random_skycoords(size=2000)
skycoords3 = get_random_skycoords(size=50000)


@pytest.mark.parametrize("skycoords", [
    skycoords1,
    skycoords2,
    skycoords3
])
def test_mocpy_from_coo_list(skycoords):
    moc = MOC.from_coo_list(skycoords, max_norder=14)


def test_moc_from_fits_image():
    import tempfile
    from astropy.io import fits
    tmp_file = tempfile.NamedTemporaryFile()
    image_path = 'notebooks/demo-data/image_with_mask.fits.gz'

    with fits.open(image_path) as hdulist:
        moc = MOC.from_image(header=hdulist[0].header,
                             moc_order=10,
                             mask_arr=hdulist[0].data)

    moc.write(tmp_file.name)


def test_simple_test_t_moc():
    t_moc = TimeMoc()
    s_time = Time(2 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    e_time = Time(7 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    t_moc.add_time_interval(s_time, e_time)
    assert t_moc.total_duration.sec == 6 * 1e-6
    assert t_moc.max_order == 29

    t_moc.write('tmoc.txt', format='json')


def test_max_order_t_moc():
    t_moc = TimeMoc()
    s_time = Time(0 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    e_time = Time(1 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    t_moc.add_time_interval(s_time, e_time)

    s_time2 = Time(3 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    e_time2 = Time(10 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    t_moc.add_time_interval(s_time2, e_time2)

    t_moc.write('tmoc2.txt', format='json')
    assert t_moc.max_order == 29


def test_from_json():
    import tempfile
    from astropy.io import fits

    tmp_file = tempfile.NamedTemporaryFile()
    image_path = 'notebooks/demo-data/image_with_mask.fits.gz'

    with fits.open(image_path) as hdulist:
        moc = MOC.from_image(header=hdulist[0].header,
                             moc_order=10)

    moc.write(tmp_file.name, format='json', write_to_file=True)

    with open(tmp_file.name, 'r') as moc_file:
        import json
        moc_d = json.load(moc_file)
        moc2 = MOC.from_json(json_moc=moc_d)
        assert moc == moc2


def test_complement_mocs():
    moc = MOC.from_moc_fits_file('notebooks/demo-data/P-GALEXGR6-AIS-FUV.fits')
    
    assert moc.complement().complement() == moc


def test_tmoc_construction():
    """
    Assert a correct t_moc loaded from a fits file is equal to the t_moc built from a CSV file
    containing a list of time intervals
    """
    time_moc = TimeMoc.from_moc_fits_file('notebooks/demo-data/TMOC/HST_SDSSg/TMoc.fits')
    time_moc2 = TimeMoc.from_csv_file(path='notebooks/demo-data/TMOC/HST_SDSSg/uniq-times.csv',
                                      delta_t=TimeMoc.order_to_time_resolution(29),
                                      format='mjd',
                                      scale='tdb')

    assert time_moc == time_moc2, 'bad tmoc construction'
