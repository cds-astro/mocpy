import pytest
from ..interval_set import IntervalSet

@pytest.fixture()
def isets():
    a = IntervalSet(intervals=[(49, 73), (53, 54), (33, 63), (65, 80), (51, 80), (100, 126), (38, 68), (61, 72), (74, 102), (27, 43)])
    b = IntervalSet(intervals=[(17, 26), (17, 41), (12, 31), (32, 61), (68, 90), (77, 105), (18, 27), (12, 35), (9, 37), (87, 97)])
    return dict(a=a, b=b)

def test_init():
    s = IntervalSet({(1, 2), (3, 5), (4, 6)})


def test_interval_set(isets):
    assert isets['a'].intervals == [(27, 126)]

    assert isets['b'].intervals == [(9, 61), (68, 105)]

    assert isets['a'].union(isets['b']).intervals == [(9, 126)]

    assert isets['a'].difference(isets['b']).intervals == [(61, 68), (105, 126)]

    assert isets['b'].difference(isets['a']).intervals == [(9, 27)]

    assert isets['a'].intersection(isets['b']).intervals == [(27, 61), (68, 105)]

    assert IntervalSet.flatten(isets['a'].intervals) == [27, 126]
