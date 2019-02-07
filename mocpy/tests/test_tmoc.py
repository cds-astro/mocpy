from astropy.time import Time
from astropy.io import ascii
from ..tmoc import TimeMOC


def test_simple_test_t_moc():
    times = Time([2/TimeMOC.DAY_MICRO_SEC, 7/TimeMOC.DAY_MICRO_SEC], format='jd', scale='tdb')
    t_moc = TimeMOC.from_times(times, delta_t=TimeMOC.order_to_time_resolution(29))
    assert t_moc.total_duration.sec == 2 * 1e-6
    assert t_moc.max_order == 29

    t_moc.write('tmoc.txt', format='json')

"""
def test_max_order_t_moc():
    t_moc = TimeMoc()
    s_time = Time(0 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    e_time = Time(1 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    t_moc.add_time_interval(s_time, e_time)

    s_time2 = Time(3 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    e_time2 = Time(10 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    t_moc.add_time_interval(s_time2, e_time2)

    t_moc.write('tmoc2.txt', format='json')
    assert t_moc.max_order == 29
"""

def test_tmoc_construction():
    """
    Assert a correct t_moc loaded from a fits file is equal to the t_moc built from a CSV file
    containing a list of time intervals
    """
    time_moc = TimeMOC.from_fits('resources/TMOC/HST_SDSSg/TMoc.fits')

    # Load HST_SDSSg from a CSV file
    data = ascii.read('resources/TMOC/HST_SDSSg/uniq-times.csv', format='csv')
    time_moc2 = TimeMOC.from_time_ranges(Time(data['t_min'], format="mjd", scale="tdb"),
                                         Time(data['t_max'], format="mjd", scale="tdb"),
                                         delta_t=TimeMOC.order_to_time_resolution(29))
    """
    time_moc2 = TimeMoc.from_csv_file(path='resources/TMOC/HST_SDSSg/uniq-times.csv',
                                      delta_t=TimeMoc.order_to_time_resolution(29),
                                      format='mjd',
                                      scale='tdb')
    """
    assert time_moc == time_moc2, 'bad tmoc construction'
