import pytest
import random
from .interval_set import IntervalSet

@pytest.fixture()
def isets():
    random.seed(0)

    a = IntervalSet()
    for x in range(0, 10):
        start = random.randint(0, 100)
        a.add((start, start + random.randint(0, 30)))

    b = IntervalSet()
    for x in range(0, 10):
        start = random.randint(0, 100)
        b.add((start, start + random.randint(0, 30)))

    return dict(a=a, b=b)


def test_interval_set(isets):

    assert isets['a'].intervals == [(27, 126)]

    assert isets['b'].intervals == [(9, 61), (68, 105)]

    assert isets['a'].union(isets['b']).intervals == [(9, 126)]

    assert isets['a'].difference(isets['b']).intervals == [(61, 68), (105, 126)]

    assert isets['b'].difference(isets['a']).intervals == [(9, 27)]

    assert isets['a'].intersection(isets['b']).intervals == [(27, 61), (68, 105)]

    assert IntervalSet.flatten(isets['a'].intervals) == [27, 126]
