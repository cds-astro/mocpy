import numpy as np

LogTable256 = np.array([
    -1,
    0,
    1, 1,
    2, 2, 2, 2,
    3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
])


def log2_lut(v):
    """
    See `this algo <https://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup>`__ for
    computing the log2 of a 32 bit integer using a look up table

    Parameters
    ----------
    v : int
        32 bit integer

    Returns
    -------

    """
    res = np.zeros(v.shape, dtype=np.int32)

    tt = v >> 16
    tt_zero = (tt == 0)
    tt_not_zero = ~tt_zero

    t_h = tt >> 8
    t_zero_h = (t_h == 0) & tt_not_zero
    t_not_zero_h = ~t_zero_h & tt_not_zero

    res[t_zero_h] = LogTable256[tt[t_zero_h]] + 16
    res[t_not_zero_h] = LogTable256[t_h[t_not_zero_h]] + 24

    t_l = v >> 8
    t_zero_l = (t_l == 0) & tt_zero
    t_not_zero_l = ~t_zero_l & tt_zero

    res[t_zero_l] = LogTable256[v[t_zero_l]]
    res[t_not_zero_l] = LogTable256[t_l[t_not_zero_l]] + 8

    return res


def uniq2orderipix_lut(uniq):
    """
    ~30% faster than the method below
    Parameters
    ----------
    uniq

    Returns
    -------

    """
    order = log2_lut(uniq >> 2) >> 1
    ipix = uniq - (1 << (2 * (order + 1)))
    return order, ipix


def uniq2orderipix(uniq):
    """
    convert a HEALPix pixel coded as a NUNIQ number
    to a (norder, ipix) tuple
    """
    order = (np.log2(uniq // np.uint8(4))) // np.uint8(2)
    order = order.astype(np.uint8)
    ipix = uniq - np.uint64(4) * (np.uint64(4) ** np.uint64(order))

    return order, ipix.astype(np.uint64)


def orderipix2uniq(n_order, n_pix):
    return n_pix + ((np.uint64(4) ** n_order) << np.uint64(2))


# x : int 64 bits
def number_trailing_zeros(x):
    x = int(x)
    bits = 0
    # convention for x == 0 => return 0
    if x == 0:
        return 0

    if not (x & 0xFFFFFFFF):
        bits += 32
        x >>= 32
    if not (x & 0xFFFF):
        bits += 16
        x >>= 16
    if not (x & 0xFF):
        bits += 8
        x >>= 8
    if not (x & 0xF):
        bits += 4
        x >>= 4
    if not (x & 0x3):
        bits += 2
        x >>= 2
    if not (x & 1):
        bits += 1
        x >>= 1

    return np.uint8(bits)
