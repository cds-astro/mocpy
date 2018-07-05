import pytest
import numpy as np
from ..interval_set import IntervalSet


@pytest.fixture()
def isets():
    a = IntervalSet(np.asarray([(49, 73), (53, 54), (33, 63), (65, 80), (51, 80), (100, 126), (38, 68), (61, 72), (74, 102), (27, 43)]))
    b = IntervalSet(np.asarray([(17, 26), (17, 41), (12, 31), (32, 61), (68, 90), (77, 105), (18, 27), (12, 35), (9, 37), (87, 97)]))
    return dict(a=a, b=b)


def test_interval_set(isets):
    assert isets['a'] == IntervalSet(np.asarray([(27, 126)]))

    assert isets['b'] == IntervalSet(np.asarray([(9, 61), (68, 105)]))

    assert isets['a'].union(isets['b']) == IntervalSet(np.asarray([(9, 126)]))

    assert isets['a'].difference(isets['b']) == IntervalSet(np.asarray([(61, 68), (105, 126)]))

    assert isets['b'].difference(isets['a']) == IntervalSet(np.asarray([(9, 27)]))

    assert isets['a'].intersection(isets['b']) == IntervalSet(np.asarray([(27, 61), (68, 105)]))
