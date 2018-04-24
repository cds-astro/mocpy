#!/usr/bin/env python3

from mocpy import MOC, TimeMoc
from astroquery.vizier import Vizier

from optparse import OptionParser
from os import path
from os import remove


def main():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-c", "--catalog",
                      action="store",
                      dest="catalog_index",
                      type=str,
                      default="I/293/npm2cros",
                      help="Name of the catalog to filter",)
    parser.add_option("-m", "--type",
                      action="store",
                      dest='type_moc',
                      type='choice',
                      default='moc',
                      choices=['moc', 'time_moc'],
                      help="Specify the type of moc")
    parser.add_option("-n", "--time_column_name",
                      action="store",
                      dest='t_column',
                      type=str,
                      default=None,
                      help="The name of the time column ")
    parser.add_option("-q", "--time_format",
                      action="store",
                      dest='time_format',
                      type=str,
                      default=None,
                      help="The time format used by the table")
    parser.add_option("-t", "--table_index",
                      action="store",
                      dest="table_index",
                      type=int,
                      default=1,
                      help="Table index of the catalog",)
    parser.add_option("-f", "--moc",
                      action="store",
                      type=str,
                      default="notebooks/demo-data/P-SDSS9-r.fits",
                      dest="filename",
                      help="Path to the fits file describing the moc",)
    parser.add_option("-p", "--plot",
                      action="store",
                      type=int,
                      default=True,
                      dest="plot",
                      help="Plot the moc built from the filtered table",)
    parser.add_option("-o", "--output_file",
                      action="store",
                      type=str,
                      default='filtered_table.fits',
                      dest="output_file",
                      help="Path to the fits output file containing the filtering table",)
    (options, args) = parser.parse_args()

    viz = Vizier(columns=['*', '_RAJ2000', '_DEJ2000'])
    viz.ROW_LIMIT = -1
    table = viz.get_catalogs(options.catalog_index)[options.table_index - 1]

    if options.type_moc is 'moc':
        moc = MOC.from_file(options.filename)
        filtered_table = moc.filter_table(table=table,
                                          ra_column='_RAJ2000',
                                          dec_column='_DEJ2000')
        result_moc = MOC.from_table(filtered_table,
                                    ra_column='_RAJ2000',
                                    dec_column='_DEJ2000',
                                    moc_order=6)
    else:
        moc = TimeMoc.from_file(options.filename)
        filtered_table = moc.filter_table(table=table,
                                          t_column=options.t_column,
                                          format=options.time_format)
        result_moc = TimeMoc.from_table(filtered_table,
                                        t_column=options.t_column,
                                        format=options.time_format)

    print('Filtered table :\n{0}'.format(filtered_table))
    if options.plot:
        result_moc.plot(title='filtered moc')

    if path.isfile(options.output_file):
        remove(options.output_file)

    filtered_table.write(options.output_file, format='votable')


if __name__ == '__main__':
    main()
