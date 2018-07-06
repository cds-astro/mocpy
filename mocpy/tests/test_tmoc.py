from astropy.time import Time
from ..tmoc import TimeMoc


def test_simple_test_t_moc():
    t_moc = TimeMoc()
    s_time = Time(2 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    e_time = Time(7 / TimeMoc.DAY_MICRO_SEC, format='jd', scale='tai')
    t_moc.add_time_interval(s_time, e_time)
    assert t_moc.total_duration.sec == 6 * 1e-6
    assert t_moc.max_order == 29

    t_moc.write('tmoc.txt', format='json')


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


def test_tmoc_construction():
    """
    Assert a correct t_moc loaded from a fits file is equal to the t_moc built from a CSV file
    containing a list of time intervals
    """
    time_moc = TimeMoc.from_fits('notebooks/demo-data/TMOC/HST_SDSSg/TMoc.fits')
    time_moc2 = TimeMoc.from_csv_file(path='notebooks/demo-data/TMOC/HST_SDSSg/uniq-times.csv',
                                      delta_t=TimeMoc.order_to_time_resolution(29),
                                      format='mjd',
                                      scale='tdb')

    assert time_moc == time_moc2, 'bad tmoc construction'
