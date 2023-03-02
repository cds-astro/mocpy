import pytest
import copy

import numpy as np

from astropy.coordinates import SkyCoord, Angle
from astropy.io.votable import parse_single_table
import astropy.units as u
from astropy.io import fits

import cdshealpix

from ..moc import MOC, WCS


@pytest.fixture()
def isets():
    a = MOC.from_depth29_ranges(
        29,
        np.array(
            [
                [49, 73],
                [53, 54],
                [33, 63],
                [65, 80],
                [51, 80],
                [100, 126],
                [38, 68],
                [61, 72],
                [74, 102],
                [27, 43],
            ],
            dtype=np.uint64,
        ),
    )
    b = MOC.from_depth29_ranges(
        29,
        np.array(
            [
                [17, 26],
                [17, 41],
                [12, 31],
                [32, 61],
                [68, 90],
                [77, 105],
                [18, 27],
                [12, 35],
                [9, 37],
                [87, 97],
            ],
            dtype=np.uint64,
        ),
    )
    return {"a": a, "b": b}


def test_interval_set_consistency(isets):
    assert isets["a"] == MOC.from_depth29_ranges(
        29,
        np.array([[27, 126]], dtype=np.uint64),
    )
    assert isets["b"] == MOC.from_depth29_ranges(
        29,
        np.array([[9, 61], [68, 105]], dtype=np.uint64),
    )


def test_uniq_hpx():
    moc = MOC.from_depth29_ranges(29, np.array([[0, 1]], dtype=np.uint64))
    uniq_hpx = np.array([4 * 4**29], dtype=np.uint64)
    assert moc.uniq_hpx == uniq_hpx

    moc = MOC.from_depth29_ranges(29, np.array([[7, 76]], dtype=np.uint64))
    uniq_hpx = np.array(
        [
            1 + 4 * 4**27,
            2 + 4 * 4**27,
            3 + 4 * 4**27,
            2 + 4 * 4**28,
            3 + 4 * 4**28,
            16 + 4 * 4**28,
            17 + 4 * 4**28,
            18 + 4 * 4**28,
            7 + 4 * 4**29,
        ],
        dtype=np.uint64,
    )
    assert (np.sort(moc.uniq_hpx) == uniq_hpx).all()


def test_to_depth29_ranges(isets):
    l = isets["a"].to_depth29_ranges
    r = np.array([[27, 126]], dtype=np.uint64)
    assert np.array_equal(l, r)
    l = isets["b"].to_depth29_ranges
    r = np.array([[9, 61], [68, 105]], dtype=np.uint64)
    assert np.array_equal(l, r)


def test_interval_set_union(isets):
    assert isets["a"].union(isets["b"]) == MOC.from_depth29_ranges(
        29,
        np.array([[9, 126]], dtype=np.uint64),
    )
    assert isets["a"].union(MOC.new_empty(29)) == MOC.from_depth29_ranges(
        29,
        np.array([[27, 126]], dtype=np.uint64),
    )
    assert MOC.new_empty(29).union(isets["a"]) == MOC.from_depth29_ranges(
        29,
        np.array([[27, 126]], dtype=np.uint64),
    )


def test_interval_set_intersection(isets):
    assert isets["a"].intersection(isets["b"]) == MOC.from_depth29_ranges(
        29,
        np.array([[27, 61], [68, 105]], dtype=np.uint64),
    )
    assert isets["a"].intersection(MOC.new_empty(29)) == MOC.new_empty(29)
    assert MOC.new_empty(29).intersection(isets["a"]) == MOC.new_empty(29)


def test_interval_set_difference(isets):
    assert isets["a"].difference(isets["b"]) == MOC.from_depth29_ranges(
        29,
        np.array([[61, 68], [105, 126]], dtype=np.uint64),
    )
    assert isets["b"].difference(isets["a"]) == MOC.from_depth29_ranges(
        29,
        np.array([[9, 27]], dtype=np.uint64),
    )
    assert MOC.new_empty(29).difference(isets["a"]) == MOC.new_empty(29)
    assert isets["a"].difference(MOC.new_empty(29)) == isets["a"]


def test_interval_set_min(isets):
    assert isets["a"].min_index == np.uint64(27)
    assert isets["b"].min_index == np.uint64(9)
    assert isets["a"].union(isets["b"]).min_index == np.uint64(9)


def test_interval_set_max(isets):
    assert isets["a"].max_index == np.uint64(126)
    assert isets["b"].max_index == np.uint64(105)
    assert isets["a"].union(isets["b"]).max_index == np.uint64(126)


def test_interval_min_depth():
    big_cells = np.array([[0, 4**29]], dtype=np.uint64)
    itv_result = MOC.from_depth29_ranges(29, big_cells)

    small_cells = np.array(
        [
            [0, 4**28],
            [4**28, 2 * 4**28],
            [2 * 4**28, 3 * 4**28],
            [3 * 4**28, 4**29],
        ],
        dtype=np.uint64,
    )
    itv_small_cells = MOC.from_depth29_ranges(29, small_cells)
    assert itv_result == itv_small_cells


def test_complement():
    assert MOC.from_depth29_ranges(
        max_depth=29,
        ranges=None,
    ).complement() == MOC.from_depth29_ranges(
        max_depth=29,
        ranges=np.array([[0, 12 * 4**29]], dtype=np.uint64),
    )
    assert MOC.new_empty(
        max_depth=29,
    ).complement().complement() == MOC.from_depth29_ranges(max_depth=29, ranges=None)
    assert MOC.from_depth29_ranges(
        29,
        np.array([[1, 2], [6, 8], [5, 6]], dtype=np.uint64),
    ).complement() == MOC.from_depth29_ranges(
        29,
        np.array([[0, 1], [2, 5], [8, 12 * 4**29]], dtype=np.uint64),
    )


# --- TESTING MOC creation ---#


def get_random_skycoords(size):
    return SkyCoord(
        ra=np.random.uniform(0, 360, size),
        dec=np.random.uniform(-90, 90, size),
        unit="deg",
    )


skycoords1 = get_random_skycoords(size=1000)
skycoords2 = get_random_skycoords(size=2000)
skycoords3 = get_random_skycoords(size=50000)


@pytest.fixture()
def skycoords_gen_f():
    def gen_f(size):
        return SkyCoord(
            np.random.uniform(0, 360, size),
            np.random.uniform(-90, 90, size),
            unit="deg",
        )

    return gen_f


@pytest.fixture()
def lonlat_gen_f():
    def gen_f(size):
        return (
            np.random.uniform(0, 360, size) * u.deg,
            np.random.uniform(-90, 90, size) * u.deg,
        )

    return gen_f


@pytest.mark.parametrize("size", [1000, 10000, 50000])
def test_moc_from_skycoords(skycoords_gen_f, size):
    skycoords = skycoords_gen_f(size)
    MOC.from_skycoords(skycoords, max_norder=7)


@pytest.mark.parametrize("size", [1000, 10000, 50000])
def test_moc_from_lonlat(lonlat_gen_f, size):
    lon, lat = lonlat_gen_f(size)
    MOC.from_lonlat(lon=lon, lat=lat, max_norder=6)


def test_from_healpix_cells():
    ipix = np.array([40, 87, 65])
    depth = np.array([3, 3, 3])
    np.array([True, True, True])

    MOC.from_healpix_cells(max_depth=29, ipix=ipix, depth=depth)


def test_moc_from_fits():
    fits_path = "resources/P-GALEXGR6-AIS-FUV.fits"
    MOC.load(fits_path, "fits")


def test_moc_from_fits_url():
    url = "http://skies.esac.esa.int/Spitzer/IRAC1_bright_ISM/Moc.fits"
    MOC.from_fits(url)


def test_moc_consistent_with_aladin():
    truth = MOC.load("resources/CDS-I-125A-catalog_MOC.fits", "fits")
    table = parse_single_table("resources/I_125A_catalog.vot").to_table()

    moc = MOC.from_lonlat(
        table["_RAJ2000"].T * u.deg,
        table["_DEJ2000"].T * u.deg,
        max_norder=8,
    )

    assert moc == truth


def test_moc_from_fits_images():
    image_path = "resources/image_with_mask.fits.gz"

    MOC.from_fits_images([image_path], max_norder=15)


def test_from_fits_images_2():
    MOC.from_fits_images(["resources/u_gal.fits"], max_norder=10)


@pytest.fixture()
def moc_from_fits_image():
    image_path = "resources/image_with_mask.fits.gz"

    with fits.open(image_path) as hdulist:
        moc = MOC.from_fits_image(hdu=hdulist[0], max_norder=7, mask=hdulist[0].data)

    return moc


@pytest.fixture()
def moc_from_json():
    return MOC.from_json({"8": [45, 78], "4": [42, 57]})


def test_moc_from_fits_image(moc_from_fits_image):
    assert isinstance(moc_from_fits_image, MOC)


def test_moc_serialize_and_from_json(moc_from_json):
    ipix_d = moc_from_json.serialize(format="json")
    moc2 = MOC.from_json(ipix_d)
    assert moc_from_json == moc2


@pytest.mark.parametrize(
    "expected, moc_str",
    [
        (
            MOC.from_json(
                {
                    "5": [8, 9, 10, 42, 43, 44, 45, 54, 46],
                    "6": [4500],
                    "7": [],
                    "8": [45],
                },
            ),
            "5/8-10 42-46 54\n\r 6/4500 8/45",
        ),
        (MOC.from_json({}), "0/"),
        (MOC.from_json({"29": [101]}), "29/101"),
        (MOC.from_json({"0": [1, 0, 9]}), "0/0-1 9"),
        (MOC.from_json({"0": [2, 9]}), "0/2 9"),
        (MOC.from_json({"0": [2], "8": [8, 9, 10], "11": []}), "0/2\r \n 8/8-10\n 11/"),
    ],
    # (MOC.from_json({"0": [2, 9], "1": [9]}), "0/2 9"),
)
def test_from_str(expected, moc_str):
    assert MOC.from_str(moc_str) == expected


@pytest.mark.parametrize(
    "expected, moc_str",
    [
        (
            MOC.from_json(
                {
                    "5": [8, 9, 10, 42, 43, 44, 45, 54, 46],
                    "6": [4500],
                    "7": [],
                    "8": [45],
                },
            ),
            "5/8-10 42-46 54\n\r 6/4500 8/45",
        ),
        (MOC.from_json({}), "0/"),
        (MOC.from_json({"29": [101]}), "29/101"),
        (MOC.from_json({"0": [1, 0, 9]}), "0/0-1 9"),
        (MOC.from_json({"0": [2, 9]}), "0/2 9"),
        (MOC.from_json({"0": [2], "8": [8, 9, 10], "11": []}), "0/2\r \n 8/8-10\n 11/"),
    ],
    # (MOC.from_json({"0": [2, 9], "1": [9]}), "0/2 9"),
)
def test_from_string(expected, moc_str):
    assert MOC.from_string(moc_str, "ascii") == expected


def test_moc_full_skyfraction():
    moc = MOC.from_json({"0": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]})
    assert moc.sky_fraction == 1.0


def test_moc_skyfraction():
    moc = MOC.from_json({"0": [0, 1, 2, 3, 4, 5]})
    assert moc.sky_fraction == 0.5


def test_sky_fraction_on_empty_coverage():
    moc = MOC.new_empty(max_depth=29)
    assert moc.sky_fraction == 0


# --- TESTING MOC serialization ---#
def test_moc_serialize_to_fits(moc_from_fits_image):
    hdulist = moc_from_fits_image.serialize(format="fits")
    assert isinstance(hdulist, fits.hdu.hdulist.HDUList)


def test_moc_serialize_to_json(moc_from_fits_image):
    moc_json = moc_from_fits_image.serialize(format="json")
    assert isinstance(moc_json, dict)


@pytest.mark.parametrize(
    "moc, expected",
    [
        (
            MOC.from_json(
                {
                    "5": [8, 9, 10, 42, 43, 44, 45, 54, 46],
                    "6": [4500],
                    "7": [],
                    "8": [45],
                },
            ),
            "5/8-10 42-46 54 6/4500 8/45 ",
        ),
        (MOC.from_json({}), "0/ "),
        (MOC.from_json({"29": [101]}), "29/101 "),
        (MOC.from_json({"0": [1, 0, 9]}), "0/0-1 9 "),
        (MOC.from_json({"0": [2, 9]}), "0/2 9 "),
    ],
    #  (MOC.from_json({"0": [2, 9], "1": [9]}), "0/2 9"),
)
def test_serialize_to_str(moc, expected):
    assert moc.serialize(format="str") == expected


@pytest.mark.parametrize(
    "filename, overwrite, format, os_error",
    [
        ("moc", True, "fits", False),
        ("moc", False, "fits", True),
        ("moc", True, "json", False),
        ("moc", True, "ascii", False),
        ("moc", False, "ascii", True),
    ],
)
def test_write(moc_from_json, filename, overwrite, format, os_error):
    if os_error:
        with pytest.raises(OSError):  # TODO add the match parameter of the exception
            moc_from_json.save(filename, format=format, overwrite=overwrite)
    else:
        moc_from_json.save(filename, format=format, overwrite=overwrite)


# --- TESTING MOC plot functions ---#
def test_mpl_fill():
    fits_path = "resources/P-GALEXGR6-AIS-FUV.fits"
    moc = MOC.load(fits_path, "fits")

    import matplotlib.pyplot as plt

    fig = plt.figure(111, figsize=(10, 10))
    with WCS(
        fig,
        fov=50 * u.deg,
        center=SkyCoord(0, 20, unit="deg", frame="icrs"),
        coordsys="icrs",
        rotation=Angle(0, u.degree),
        projection="AIT",
    ) as wcs:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc.fill(ax=ax, wcs=wcs, alpha=0.5, color="r")


def test_mpl_border():
    fits_path = "resources/P-GALEXGR6-AIS-FUV.fits"
    moc = MOC.load(fits_path, "fits")

    import matplotlib.pyplot as plt

    fig = plt.figure(111, figsize=(10, 10))
    with WCS(
        fig,
        fov=50 * u.deg,
        center=SkyCoord(0, 20, unit="deg", frame="icrs"),
        coordsys="icrs",
        rotation=Angle(0, u.degree),
        projection="AIT",
    ) as wcs:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc.border(ax=ax, wcs=wcs, color="g")


# --- TESTING MOC features ---#
@pytest.mark.parametrize("order", [4, 5, 6, 15, 20, 28])
def test_moc_contains(order):
    # defines 20 random healpix cells of the required order
    size = 20
    healpix_arr = np.random.randint(0, 12 * 4**order, size, dtype="uint64")
    # defines a moc containing the 20 points
    moc = MOC.from_json(json_moc={str(order): np.unique(healpix_arr).tolist()})
    # the complementary should not contain them
    moc_complement = moc.complement()
    # coordinates of the 20 random points
    lon, lat = cdshealpix.healpix_to_lonlat(healpix_arr, order)
    # tests
    should_be_inside_arr = moc.contains_lonlat(lon=lon, lat=lat)
    assert should_be_inside_arr.all()
    should_be_outside_arr = moc_complement.contains_lonlat(lon=lon, lat=lat)
    assert not should_be_outside_arr.any()
    # test keep_inside field
    should_be_outside_arr = moc.contains_lonlat(lon=lon, lat=lat, keep_inside=False)
    assert not should_be_outside_arr.any()
    should_be_inside_arr = moc_complement.contains_lonlat(
        lon=lon,
        lat=lat,
        keep_inside=False,
    )
    assert should_be_inside_arr.all()


# test 2d-arrays as lon lat input
def test_moc_contains_2d_parameters():
    """Test that not only 1d arrays are accepted."""
    lon = Angle(np.array([[1, 2, 3], [2, 40, -5]]), unit=u.deg)
    lat = Angle(np.array([[20, 25, 10], [-60, 80, 0]]), unit=u.deg)
    lat2 = Angle(np.array([[20, 25, 10, 22], [-60, 80, 0, 10]]), unit=u.deg)
    moc = MOC.from_polygon(lon=lon, lat=lat, max_depth=12)
    should_be_inside = moc.contains_lonlat(lon=lon, lat=lat)
    assert should_be_inside.all()
    # test mismatched
    with pytest.raises(ValueError, match=r"Lon shape different from lat shape.*"):
        moc.contains_lonlat(lon=lon, lat=lat2)


def test_degrade_to_order():
    hst_fits_path = "resources/hst.fits"
    hst_moc = MOC.load(hst_fits_path, "fits")

    max_depth = hst_moc.max_order

    for order in reversed(range(0, max_depth)):
        hst_moc = hst_moc.degrade_to_order(order)
        assert hst_moc.sky_fraction <= 1.0


def test_from_ring():
    MOC.from_ring(
        lon=0 * u.deg,
        lat=0 * u.deg,
        internal_radius=Angle(5, u.deg),
        external_radius=Angle(10, u.deg),
        max_depth=10,
    )


# TODO: IMPROVE THE ALGO
"""
def test_boundaries():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.load(fits_path, 'fits')
    moc = moc.degrade_to_order(6)
    boundaries_l = moc.get_boundaries()
"""


def test_from_elliptical_cone():
    MOC.from_elliptical_cone(
        lon=0 * u.deg,
        lat=0 * u.deg,
        a=Angle(10, u.deg),
        b=Angle(5, u.deg),
        pa=Angle(0, u.deg),
        max_depth=10,
    )


@pytest.fixture()
def mocs():
    moc1 = {"1": [0]}
    moc1_increased = {"0": [0], "1": [17, 19, 22, 23, 35]}
    moc2 = {"1": [30]}
    moc2_increased = {"0": [7], "1": [8, 9, 25, 43, 41]}

    return {
        "moc1": MOC.from_json(moc1),
        "moc1_increased": MOC.from_json(moc1_increased),
        "moc2": MOC.from_json(moc2),
        "moc2_increased": MOC.from_json(moc2_increased),
    }


def test_add_neighbours(mocs):
    print(mocs["moc1"])
    mocs["moc1"].add_neighbours()
    print(mocs["moc1"])
    assert mocs["moc1"] == mocs["moc1_increased"]

    mocs["moc2"].add_neighbours()
    assert mocs["moc2"] == mocs["moc2_increased"]


def test_remove_neighbours(mocs):
    mocs["moc1_increased"].remove_neighbours()
    mocs["moc2_increased"].remove_neighbours()
    assert mocs["moc1_increased"] == mocs["moc1"]
    assert mocs["moc2_increased"] == mocs["moc2"]


def test_neighbours(mocs):
    moc1 = copy.deepcopy(mocs["moc1"])
    moc2 = copy.deepcopy(mocs["moc2"])
    moc1.add_neighbours().remove_neighbours()
    moc2.add_neighbours().remove_neighbours()
    assert moc1 == mocs["moc1"]
    assert moc2 == mocs["moc2"]


# --- TESTING MOC operations ---#
@pytest.fixture()
def mocs_op():
    moc1 = MOC.from_json({"0": [0, 2, 3, 4, 5]})
    moc2 = MOC.from_json({"0": [0, 1, 7, 4, 3]})
    return {"first": moc1, "second": moc2}


def test_moc_union(mocs_op):
    assert mocs_op["first"].union(mocs_op["second"]) == MOC.from_json(
        {"0": [0, 1, 2, 3, 4, 5, 7]},
    )
    assert mocs_op["first"] + mocs_op["second"] == MOC.from_json(
        {"0": [0, 1, 2, 3, 4, 5, 7]},
    )
    assert mocs_op["first"] | mocs_op["second"] == MOC.from_json(
        {"0": [0, 1, 2, 3, 4, 5, 7]},
    )


def test_moc_intersection(mocs_op):
    assert mocs_op["first"].intersection(mocs_op["second"]) == MOC.from_json(
        {"0": [0, 3, 4]},
    )
    assert mocs_op["first"] & mocs_op["second"] == MOC.from_json({"0": [0, 3, 4]})


def test_moc_difference(mocs_op):
    assert mocs_op["first"].difference(mocs_op["second"]) == MOC.from_json(
        {"0": [2, 5]},
    )
    assert mocs_op["first"] - mocs_op["second"] == MOC.from_json({"0": [2, 5]})


def test_moc_complement_consistency():
    moc = MOC.load("resources/P-GALEXGR6-AIS-FUV.fits", "fits")
    assert moc.complement().complement() == moc


def test_from_fits_old():
    moc = MOC.from_fits("resources/V_147_sdss12.moc.fits")
    assert moc.complement().complement() == moc


@pytest.mark.parametrize(
    "input, expected",
    [
        (
            MOC.from_json({"0": [1, 3]}),
            MOC.from_json({"0": [0, 2, 4, 5, 6, 7, 8, 9, 10, 11]}),
        ),
    ],
)
def test_moc_complement(input, expected):
    assert input.complement() == expected
    assert ~input == expected


def test_spatial_res_to_order():
    order = np.arange(14)

    res = MOC.order_to_spatial_resolution(order)
    output = MOC.spatial_resolution_to_order(res)

    assert (order == output).all()


def test_from_valued_healpix_cells_empty():
    uniq = np.array([])
    values = np.array([])

    MOC.from_valued_healpix_cells(uniq, values, 12)


def test_from_valued_healpix_cells_different_sizes():
    uniq = np.array([500])
    values = np.array([])

    with pytest.raises(
        ValueError,
        match="`uniq` and values do not have the same size.",
    ):
        MOC.from_valued_healpix_cells(uniq, values, 12)


def test_from_valued_healpix_cells_cumul_from_sup_cumul_to():
    uniq = np.array([500])
    values = np.array([1.0])

    with pytest.raises(ValueError, match="`cumul_from` has to be < to `cumul_to`."):
        MOC.from_valued_healpix_cells(uniq, values, 12, cumul_from=0.8, cumul_to=-5.0)


@pytest.mark.parametrize(
    "cumul_from, cumul_to",
    [(-5.0, 1.0), (np.nan, np.inf), (np.nan, np.nan), (np.inf, np.nan), (-10.0, -5.0)],
)
def test_from_valued_healpix_cells_weird_values(cumul_from, cumul_to):
    uniq = np.array([500])
    values = np.array([-1.0])

    MOC.from_valued_healpix_cells(
        uniq,
        values,
        12,
        cumul_from=cumul_from,
        cumul_to=cumul_to,
    )


def test_from_valued_healpix_cells_bayestar():
    from astropy.io import fits

    fits_image_filename = "./resources/bayestar.multiorder.fits"

    with fits.open(fits_image_filename) as hdul:
        hdul.info()
        hdul[1].columns

        data = hdul[1].data

    uniq = data["UNIQ"]
    probdensity = data["PROBDENSITY"]

    import astropy_healpix as ah
    import astropy.units as u

    level, _ = ah.uniq_to_level_ipix(uniq)
    area = ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.steradian)

    prob = probdensity * area

    cumul_to = np.linspace(0.01, 2.0, num=10)

    for b in cumul_to:
        MOC.from_valued_healpix_cells(uniq, prob, 12, cumul_from=0.0, cumul_to=b)


def test_from_valued_healpix_cells_bayestar_and_split():
    fits_mom_filename = "./resources/bayestar.multiorder.fits"
    moc = MOC.from_multiordermap_fits_file(
        fits_mom_filename,
        cumul_from=0.0,
        cumul_to=0.9,
    )
    count = moc.split_count()
    assert count == 2
    mocs = list(moc.split())
    assert len(mocs) == 2
    for moc in mocs:
        assert moc.max_order == 11


# --- TESTING new features ---#
def test_moc_save_load_deser():
    smoc = MOC.from_string("3/3 10 4/16-18 22 5/19-20 17/222 28/123456789 29/", "ascii")
    smoc_ascii = smoc.to_string("ascii")
    smoc_ascii
    smoc_json = smoc.to_string("json")
    smoc_json
    smoc_bis = MOC.from_string(smoc_json, "json")
    assert smoc == smoc_bis

    smoc_bis = MOC.load("resources/MOC2.0/smoc.ascii.txt", "ascii")
    assert smoc == smoc_bis

    smoc_bis = MOC.load("resources/MOC2.0/SMOC.fits", "fits")
    assert smoc == smoc_bis

    smoc.save("resources/MOC2.0/smoc.py.test.fits", format="fits", overwrite=True)
    smoc.save("resources/MOC2.0/smoc.py.test.json", format="json", overwrite=True)
    smoc.save("resources/MOC2.0/smoc.py.test.ascii", format="ascii", overwrite=True)
    smoc_bis = MOC.load("resources/MOC2.0/smoc.py.test.fits", "fits")
    assert smoc == smoc_bis
    smoc_bis = MOC.load("resources/MOC2.0/smoc.py.test.json", "json")
    assert smoc == smoc_bis
    smoc_bis = MOC.load("resources/MOC2.0/smoc.py.test.ascii", "ascii")
    assert smoc == smoc_bis
