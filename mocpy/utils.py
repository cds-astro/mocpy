import numpy as np


def uniq2orderipix(uniq):
    """
    convert a HEALPix pixel coded as a NUNIQ number
    to a (norder, ipix) tuple
    """
    order = ((np.log2(uniq//4)) // 2)
    order = order.astype(int)
    ipix = uniq - 4 * (4**order)
    
    return order, ipix


def orderipix2uniq(n_order, n_pix):
    return n_pix + ((4**n_order) << 2)


# x : int 64 bits
def number_trailing_zeros(x):
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

    return bits

"""
def number_trailing_zeros(i):
    # O(log(i)) complexity here. Can be done in O(1) using more low-level algo
    nb = 0
    while (i & 1) == 0 and i > 0:
        nb += 1
        i = i >> 1

    return nb
"""
