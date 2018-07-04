import pytest
import numpy as np
from astropy.coordinates import SkyCoord

from ..moc import MOC


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
def test_moc_from_skycoords(skycoords):
    moc = MOC.from_skycoords(skycoords, max_norder=7)


def test_moc_from_fits():
    fits_path = 'notebooks/demo-data/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.from_fits(fits_path)


def test_moc_from_fits_image():
    import tempfile
    from astropy.io import fits
    tmp_file = tempfile.NamedTemporaryFile()
    image_path = 'notebooks/demo-data/image_with_mask.fits.gz'

    with fits.open(image_path) as hdulist:
        moc = MOC.from_image(header=hdulist[0].header,
                             max_norder=10,
                             mask_arr=hdulist[0].data)

    moc.write(tmp_file.name)


def test_moc_from_json():
    import tempfile
    from astropy.io import fits

    tmp_file = tempfile.NamedTemporaryFile()
    image_path = 'notebooks/demo-data/image_with_mask.fits.gz'

    with fits.open(image_path) as hdulist:
        moc = MOC.from_image(header=hdulist[0].header,
                             max_norder=10)

    moc.write(tmp_file.name, format='json', write_to_file=True)

    with open(tmp_file.name, 'r') as moc_file:
        import json
        moc_d = json.load(moc_file)
        moc2 = MOC.from_json(json_moc=moc_d)
        assert moc == moc2


def test_moc_complement():
    moc = MOC.from_fits('notebooks/demo-data/P-GALEXGR6-AIS-FUV.fits')

    assert moc.complement().complement() == moc
