import pytest
import numpy as np
from ..interval_set import IntervalSet


@pytest.fixture()
def isets():
    a = IntervalSet(np.array([[49, 73], [53, 54], [33, 63], [65, 80],
        [51, 80], [100, 126], [38, 68], [61, 72],
        [74, 102], [27, 43]], dtype=np.uint64))
    b = IntervalSet(np.array([[17, 26], [17, 41], [12, 31], [32, 61],
        [68, 90], [77, 105], [18, 27], [12, 35],
        [9, 37], [87, 97]], dtype=np.uint64))
    return dict(a=a, b=b)


def test_interval_set_consistency(isets):
    assert isets['a'] == IntervalSet(np.array([[27, 126]], dtype=np.uint64))
    assert isets['b'] == IntervalSet(np.array([[9, 61], [68, 105]], dtype=np.uint64))


def test_interval_min_depth():
    big_cells = np.array([[0, 4**29]], dtype=np.uint64)
    itv_result = IntervalSet(big_cells, min_depth=1)
    
    small_cells = np.array([[0, 4**28], [4**28, 2*4**28], [2*4**28, 3*4**28], [3*4**28, 4**29]], dtype=np.uint64)
    itv_small_cells = IntervalSet(small_cells, make_consistent=False)
    assert itv_result == itv_small_cells


def test_interval_set_union(isets):
    assert isets['a'].union(isets['b']) == IntervalSet(np.array([[9, 126]], dtype=np.uint64))
    assert isets['a'].union(IntervalSet()) == IntervalSet(np.array([[27, 126]], dtype=np.uint64))
    assert IntervalSet().union(isets['a']) == IntervalSet(np.array([[27, 126]], dtype=np.uint64))


def test_interval_set_intersection(isets):
    assert isets['a'].intersection(isets['b']) == IntervalSet(np.array([[27, 61], [68, 105]], dtype=np.uint64))
    assert isets['a'].intersection(IntervalSet()) == IntervalSet()
    assert IntervalSet().intersection(isets['a']) == IntervalSet()


def test_interval_set_difference(isets):
    assert isets['a'].difference(isets['b']) == IntervalSet(np.array([[61, 68], [105, 126]], dtype=np.uint64))
    assert isets['b'].difference(isets['a']) == IntervalSet(np.array([[9, 27]], dtype=np.uint64))
    assert IntervalSet().difference(isets['a']) == IntervalSet()
    assert isets['a'].difference(IntervalSet()) == isets['a']


def test_interval_set_complement():
    assert IntervalSet().complement() == IntervalSet(np.array([[0, 12*4**29]], dtype=np.uint64))
    assert IntervalSet().complement().complement() == IntervalSet()
    assert IntervalSet(np.array([[1, 2], [6, 8], [5, 6]], dtype=np.uint64)).complement() == \
        IntervalSet(np.array([[0, 1], [2, 5], [8, 12*4**29]], dtype=np.uint64))


@pytest.fixture()
def isets2():
    nested1 = IntervalSet(np.array([[0, 1]], dtype=np.uint64))
    nuniq1 = np.array([4*4**29], dtype=np.uint64)
    nested2 = IntervalSet(np.array([[7, 76]], dtype=np.uint64))
    nuniq2 = np.array([1 + 4*4**27, 2 + 4*4**27, 3 + 4*4**27,
                      2 + 4*4**28, 3 + 4*4**28,
                      16 + 4*4**28, 17 + 4*4**28, 18 + 4*4**28,
                      7 + 4*4**29], dtype=np.uint64)
    return {
        'nest1': nested1,
        'uniq1': nuniq1,
        'nest2': nested2,
        'uniq2': nuniq2,
    }


def test_to_uniq(isets2):
    assert (isets2['nest1'].uniq == isets2['uniq1']).all()
    assert (isets2['nest2'].uniq == isets2['uniq2']).all()
    # empty nested interval set
    assert (IntervalSet().uniq == np.array([], dtype=np.uint64)).all()


def test_from_uniq(isets2):
    assert IntervalSet.from_uniq(isets2['uniq1']) == isets2['nest1']
    assert IntervalSet.from_uniq(isets2['uniq2']) == isets2['nest2']
    # empty nuniq interval set
    assert IntervalSet.from_uniq(np.array([], dtype=np.uint64)) == IntervalSet()


def test_from_to_interval_set(isets2):
    assert IntervalSet.from_uniq(isets2['nest1'].uniq) == isets2['nest1']


def test_interval_set_min(isets):
    assert isets['a'].min == np.uint64(27)
    assert isets['b'].min == np.uint64(9)
    assert isets['a'].union(isets['b']).min == np.uint64(9)


def test_interval_set_max(isets):
    assert isets['a'].max == np.uint64(126)
    assert isets['b'].max == np.uint64(105)
    assert isets['a'].union(isets['b']).max == np.uint64(126)


def test_repr_interval_set(isets):
    assert repr(isets['a']) == "[[ 27 126]]"
    assert repr(isets['b']) == "[[  9  61]\n" \
                               " [ 68 105]]"