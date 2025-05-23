import re
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest

from ..fmoc import FrequencyMOC as FMOC
from ..moc import MOC
from ..sfmoc import SFMOC


def test_new_empty_sfmoc():
    sfmoc = SFMOC.new_empty(20, 12)
    assert sfmoc.is_empty()
    assert sfmoc.max_order == (20, 12)

    # too big order
    assert SFMOC.new_empty(0, 29) == SFMOC.new_empty(0, 100)


def test_equality_on_empty():
    sfmoc = SFMOC.new_empty(20, 12)
    sfmoc2 = SFMOC.new_empty(20, 12)
    assert sfmoc == sfmoc2


def test_n_cells():
    assert SFMOC.n_cells(0, dimension="space") == 12
    i = 15
    assert SFMOC.n_cells(i, dimension="frequency") == 2 * 2**i
    # test the errors
    with pytest.raises(
        ValueError, match="The depth should be comprised between 0 and 29, but*"
    ):
        SFMOC.n_cells(-1, dimension="space")
    with pytest.raises(ValueError, match="Dimension should be either *"):
        SFMOC.n_cells(0, dimension="blah")


def test_shape_mismatch():
    freqs = [1, 2] * u.Hz
    lon = [0, 1, 2] * u.deg
    lat = [0, 1, 2] * u.deg
    with pytest.raises(
        ValueError, match="Frequencies and positions must have the same lengths."
    ):
        SFMOC.from_frequencies_and_positions(
            freqs, lon, lat, max_order_frequency=20, max_order_space=12
        )


def test_from_frequencies_and_positions():
    frequencies = [1, 2, 3] * u.Hz
    lon = [0, 1, 2] * u.deg
    lat = [0, 1, 2] * u.deg
    sfmoc = SFMOC.from_frequencies_and_positions(
        frequencies, lon, lat, max_order_frequency=20, max_order_space=12
    )
    assert not sfmoc.is_empty()
    assert sfmoc.min_frequency <= 1 * u.Hz
    assert sfmoc.max_frequency >= 3 * u.Hz

    # scalar case
    sfmoc = SFMOC.from_frequencies_and_positions(
        0.5 * u.Hz, 0 * u.deg, 0 * u.deg, max_order_frequency=0, max_order_space=0
    )
    assert isinstance(sfmoc, SFMOC)

    # errors
    with pytest.raises(
        ValueError,
        match="Frequencies and positions must be scalar quantities or 1D arrays.",
    ):
        SFMOC.from_frequencies_and_positions(
            [[0.01], [0.01]] * u.Hz,
            [[1], [1]] * u.deg,
            [[1], [1]] * u.deg,
            max_order_frequency=0,
            max_order_space=0,
        )


def test_from_frequency_ranges_and_positions():
    fmin = [0.01, 0.02, 0.03] * u.Hz
    fmax = [0.1, 0.2, 0.3] * u.Hz
    lon = [0, 1, 2] * u.deg
    lat = [0, 1, 2] * u.deg
    sfmoc = SFMOC.from_frequency_ranges_and_positions(
        fmin, fmax, lon, lat, max_order_frequency=5, max_order_space=5
    )
    assert all(sfmoc.contains(fmin, lon, lat))


def test_query_by_frequency():
    sfmoc = SFMOC.from_string("f10/0-20\ns12/0-100", format="ascii")
    fmoc = FMOC.from_string("10/0-5")
    result = sfmoc.query_by_frequency(fmoc)
    assert result == MOC.from_string("12/0-100")


def test_query_by_space():
    sfmoc = SFMOC.from_string("f10/0-20\ns12/0-100", format="ascii")
    smoc = MOC.from_string("12/0-100")
    result = sfmoc.query_by_space(smoc)
    assert result == FMOC.from_string("10/0-20")


def test_contains():
    moc = MOC.from_cone(0 * u.deg, 0 * u.deg, radius=10 * u.deg, max_depth=10)
    sfmoc = SFMOC.from_spatial_coverages(
        0.01 * u.Hz, 100 * u.Hz, moc, max_order_frequency=20
    )
    result = sfmoc.contains([0.1, 1000] * u.Hz, [0, 5] * u.deg, [0, 5] * u.deg)
    assert all(result == np.array([True, False]))


@pytest.mark.parametrize(
    ("file_name", "format_name"),
    [("sfmoc.fits", "fits"), ("sfmoc.json", "json"), ("sfmoc.ascii", "ascii")],
)
def test_roundtrip_io(file_name, format_name, tmp_path):
    sfmoc = SFMOC.from_string("f10/0-10\ns10/0-10", format="ascii")
    path = tmp_path / file_name
    sfmoc.save(path, format=format_name)
    loaded = SFMOC.load(path, format=format_name)
    assert loaded == sfmoc


def test_wrong_format_io(tmp_path):
    sfmoc = SFMOC.from_string("f10/0-10\ns10/0-10", format="ascii")
    path = tmp_path / "test_wrong_format.fits"
    sfmoc.save(path)
    with pytest.raises(ValueError, match="format should be *"):
        SFMOC.load(path, format="blah")


def test_read_fits_file():
    sfmoc = SFMOC.load(
        Path(__file__).parent.parent.parent.parent
        / "resources"
        / "SFMOC"
        / "SFMocMUSE.fits"
    )
    assert isinstance(sfmoc, SFMOC)
    assert not sfmoc.is_empty()


def test_from_string():
    sfmoc_json = SFMOC.from_string(
        """[{"f": {
        "7": [0],
        "9": [4],
        "10": [10]
    },
    "s": {
        "9": [0, 1],
        "10": [8, 9, 10]
    }
    },
    { "f": { "10": [] }, "s": { "10": [] } }
    ]
    """,
        format="json",
    )
    sfmoc_string = SFMOC.from_string("f10/0-10\ns10/0-10", format="ascii")
    assert sfmoc_json == sfmoc_string

    # error on strange format
    with pytest.raises(
        ValueError, match=re.escape("format should be one of ('ascii', 'json')")
    ):
        SFMOC.from_string("f10/0-10\ns10/0-10", format="blah")
