import numpy as np

def uniq2orderipix(uniq):
    """
    convert a HEALPix pixel coded as a NUNIQ number
    to a (norder, ipix) tuple
    """
    order = (np.log2(uniq // np.uint8(4))) // np.uint8(2)
    order = order.astype(np.uint8)
    ipix = uniq - np.uint64(4) * (np.uint64(4) ** np.uint64(order))

    return order, ipix.astype(np.uint64)
