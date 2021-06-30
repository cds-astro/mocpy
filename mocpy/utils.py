import numpy as np
from astropy.time import Time

DAY_MICRO_SEC = 86400000000.

def uniq2orderipix(uniq):
    """
    convert a HEALPix pixel coded as a NUNIQ number
    to a (norder, ipix) tuple
    """
    order = (np.log2(uniq // np.uint8(4))) // np.uint8(2)
    order = order.astype(np.uint8)
    ipix = uniq - np.uint64(4) * (np.uint64(4) ** np.uint64(order))

    return order, ipix.astype(np.uint64)

def times_to_microseconds(times):
    """
    Convert a `astropy.time.Time` into an array of integer microseconds since JD=0, keeping
    the microsecond resolution required for `~mocpy.tmoc.TimeMOC`.

    Parameters
    ----------
    times : `astropy.time.Time`
    Astropy observation times

    Returns
    -------
    times_microseconds : `np.array`
    """
    times_jd = np.asarray(times.jd, dtype=np.uint64)
    times_us = np.asarray((times - Time(times_jd, format='jd', scale='tdb')).jd * DAY_MICRO_SEC, dtype=np.uint64)
    
    return times_jd * np.uint64(DAY_MICRO_SEC) + times_us

def microseconds_to_times(times_microseconds):
    """
    Convert an array of integer microseconds since JD=0, to an array of `astropy.time.Time`.

    Parameters
    ----------
    times_microseconds : `np.array`

    Returns
    -------
    times : `astropy.time.Time`
    """
    jd1 = np.asarray(times_microseconds // DAY_MICRO_SEC, dtype=np.float64)
    jd2 = np.asarray((times_microseconds - jd1 * DAY_MICRO_SEC) / DAY_MICRO_SEC, dtype=np.float64)

    return Time(val=jd1, val2=jd2, format='jd', scale='tdb')
