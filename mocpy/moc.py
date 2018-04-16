#!/usr/bin/env python
# -*- coding: utf-8 -*

"""

moc.py:
  functions to read/write and manipulate MOCs

"""

__author__ = "Thomas Boch"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

try:
    set
except NameError:
    from sets import Set as set
import numpy as np
from .interval_set import IntervalSet
from . import utils

from astropy.table import Table
from astropy.coordinates import SkyCoord, ICRS
from astropy_healpix.healpy import ang2pix
from astropy_healpix.healpy import ring2nest
from astropy.io import fits
from astropy import wcs
from astropy_healpix import HEALPix
from astropy import units as u

import sys
if sys.version > '3':
    long = int

# Python 3 support
try:
    xrange
except NameError:
    xrange = range


def number_trailing_zeros(i):
    # TODO : there must be a smarter way of doing that
    nb = 0
    while (i & 1) == 0:
        nb += 1
        i >>= 1

    return nb


class MOC:
    HPY_MAX_NORDER = 29 # upper norder limit for healpy
    VIZ_TABLE_MOC_ROOT_URL = ''
    VIZ_CAT_MOC_ROOT_URL = ''
    
    def __init__(self):
        self._interval_set = IntervalSet() # set of intervals at HPY_MAX_NORDER norder
        self._counter = 0
        self._order = None

    def __repr__(self):
        return self._interval_set.__repr__()

    def __eq__(self, another_moc):
        if not isinstance(another_moc, MOC):
            raise TypeError
        return self._interval_set == another_moc._interval_set and\
               self._order == another_moc._order

    @property
    def max_order(self):
        """
        This returns the deepest order needed to describe the current _interval_set
        """
        # TODO: cache value
        combo = long(0)
        for iv in self._interval_set.intervals:
            combo |= iv[0] | iv[1]

        ret = MOC.HPY_MAX_NORDER - int(number_trailing_zeros(combo)//2)
        if ret < 0:
            ret = 0

        return ret

    def intersection(self, another_moc):
        """
        intersection with another MOC
        """
        iv_set_intersection = self._interval_set.intersection(another_moc._interval_set)
        
        return MOC.from_interval_set(iv_set_intersection)
    
    def union(self, another_moc):
        """
        union with another MOC
        """
        iv_set_union = self._interval_set.union(another_moc._interval_set)
        
        return MOC.from_interval_set(iv_set_union)

    def degrade_to_order(self, new_order):    
        shift = 2 * (MOC.HPY_MAX_NORDER - new_order)
        ofs = (long(1) << shift) - 1
        mask = ~ofs
        adda = long(0)
        addb = ofs
        iv_set = IntervalSet()
        for iv in self._interval_set.intervals:
            a = (iv[0] + adda) & mask
            b = (iv[1] + addb) & mask
            if b > a:
                iv_set.add((a, b))
               
        m = MOC.from_interval_set(iv_set)
        m._order = new_order
        return m

    @property
    def sky_fraction(self):
        """
        return the sky fraction (between 0 and 1) covered by the MOC
        """
        nb_pix_filled = 0
        for val in self.best_res_pixels_iterator():
            nb_pix_filled += 1

        return nb_pix_filled / float(12 * 4**self.max_order)

    def plot(self, title='MOC', coord='C'):
        """
        plot current instance using matplotlib
        coord can be 'C', 'G' or 'E' (respectively celestial, galactic or ecliptic)
        """
        from matplotlib import pyplot as plt
        import healpy as hp
        # degrade to NORDER 8 if norder is greater
        if self.max_order>8:
            plotted_moc = self.degrade_to_order(8)
        else:
            plotted_moc = self
            
        m = np.zeros(hp.nside2npix(2 ** plotted_moc.max_order))
        for val in plotted_moc.best_res_pixels_iterator():
            m[val] = 1
        
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list('w2r', ['#ffffff', '#ff0000'])
        cmap.set_under('w') 
        cmap.set_bad('gray')
        hp.mollview(m, nest=True, coord=['C', coord], title=title, cbar=False, cmap=cmap)
        hp.graticule()
        plt.show()

    def query_simbad(self, max_rows=10000):
        """
        query a view of SIMBAD data
        for SIMBAD objects in the coverage of the MOC instance
        """
        return self._query('SIMBAD', max_rows)

    def query_vizier_table(self, table_id, max_rows=10000):
        """
        query a VizieR table
        for sources in the coverage of the MOC instance
        """
        return self._query(table_id, max_rows)

    def _query(self, resource_id, max_rows):
        """
        internal method to query Simbad or a VizieR table
        for sources in the coverage of the MOC instance
        """
        import requests
        import tempfile, os

        if max_rows is not None and max_rows>=0:
            max_rows_str = str(max_rows)
        else:
            max_rows_str = str(9999999999)

        tmp_moc = tempfile.NamedTemporaryFile(delete = False)
        self.write(tmp_moc.name)
        r = requests.post('http://cdsxmatch.u-strasbg.fr/QueryCat/QueryCat', data={'mode': 'mocfile' , 'catName': resource_id, 'format': 'votable', 'limit': max_rows_str}, files={'moc': open(tmp_moc.name, 'rb')}, headers={'User-Agent': 'MOCPy'}, stream=True)
        
        tmp_vot = tempfile.NamedTemporaryFile(delete = False)
        with open(tmp_vot.name, 'w') as h:
            for line in r.iter_lines():
                if line:
                    h.write(line.decode(r.encoding)+'\n')

        from astropy.io.votable import parse_single_table
        table = parse_single_table(tmp_vot.name).to_table()

        # finally delete temp files
        os.unlink(tmp_moc.name)
        os.unlink(tmp_vot.name)

        return table

    @classmethod
    def from_file(cls, local_path, moc_order=None):
        return MOC_io.read_local(path=local_path, moc_order=moc_order)

    @classmethod
    def from_vizier_table(cls, table_id, nside=256):
        """
        return the MOC for a given VizieR table
        """
        if nside not in (8, 16, 32, 64, 128, 256, 512):
            raise Exception('Bad value for nside')
        
        return cls.from_ivorn('ivo://CDS/' + table_id, nside)
    
    MOC_SERVER_ROOT_URL = 'http://alasky.unistra.fr/MocServer/query'        
        
    @classmethod
    def from_ivorn(cls, ivorn, nside=256):
        """
        return the MOC for a given IVORN
        """
        try:
            from urllib import urlencode
        except ImportError:
            from urllib.parse import urlencode

        return cls.from_url('%s?%s' % (MOC.MOC_SERVER_ROOT_URL, urlencode({'ivorn': ivorn, 'get': 'moc', 'order': int(np.log2(nside))})))

    @classmethod        
    def from_url(cls, url):
        from astropy.utils.data import download_file

        path = download_file(url, show_progress=False, timeout=60)
        return cls.from_file(path)
        
    def add_uniq_ipix(self, uniq_ipix):
        """
        add a UNIQ pixel to the current MOC
        """
        order, ipix = utils.uniq2orderipix(uniq_ipix)
        self.add_pix(order, ipix)
       
    def add_pix_list(self, order, ipix_list, nest=True):
        """
        add a list of HEALPix pixel number at a given order
        to the current object
        """
        for ipix in ipix_list:
            self.add_pix(order, ipix, nest)

    def _ensure_consistency(self):
        """
        Force IntervalSet internal consistency
        """
        self._interval_set.intervals

    def add_pix(self, order, ipix, nest=True):
        """
        add a given HEALPix pixel number to the current object
        """
        self._counter += 1
        # force consistency to prevent too large interval array
        if self._counter == 1000:
            self._ensure_consistency()
            self._counter = 0

        if order > MOC.HPY_MAX_NORDER:
            raise Exception('norder can not be greater than MOC max norder')
        
        if not nest:
            ipix = ring2nest(1 << order, ipix)

        p1 = ipix
        p2 = ipix + 1
        shift = 2 * (MOC.HPY_MAX_NORDER - order)

        self._interval_set.add((p1 << shift, p2 << shift))

    def to_uniq_interval_set(self):
        """
        Return an IntervalSet in NUNIQ norder
        reflecting the current state of the MOC
        """
        r2 = IntervalSet(self._interval_set)
        r3 = IntervalSet()
        res = IntervalSet()

        intervals = self._interval_set.intervals
        if len(intervals) == 0:
            return res
        
        max_res_order = MOC.HPY_MAX_NORDER

        for order in xrange(0, max_res_order):
            if r2.empty():
                return res
            
            shift = long(2 * (max_res_order - order))
            ofs = (long(1) << shift)-1
            ofs2 = long(1) << (2*order+2)
            
            r3.clear()
            for iv in r2.intervals:
                a = (long(iv[0]) + ofs) >> shift
                b = iv[1] >> shift
                r3.add((a << shift, b << shift))
                res.add((a + ofs2, b + ofs2))
            
            if not r3.empty():
                r2 = r2.difference(r3)
      
        return res

    def add_position(self, ra, dec, max_norder):
        """
        add the HEALPix bin containing the (ra, dec) position
        """
        theta, phi = utils.radec2thetaphi(ra, dec)
        ipix = ang2pix(2**max_norder, theta, phi, nest=True)

        if isinstance(ipix, np.ndarray):
            for i in range(ipix.shape[0]):
                self.add_pix(max_norder, ipix[i])
        else:
            self.add_pix(max_norder, ipix)

        self._order = max_norder

    @classmethod
    def from_json(cls, json_moc):
        uniq_interval = IntervalSet()
        for n_order, n_pix_l in json_moc.items():
            n_order = int(n_order)

            uniq_interval_add = np.vectorize(uniq_interval.add)
            uniq_interval_add(utils.orderipix2uniq(n_order, np.array(n_pix_l)))

        return MOC.from_uniq_interval_set(uniq_interval)

    @classmethod
    def from_uniq_interval_set(cls, uniq_is):
        """
        Create a MOC from an IntervalSet of NUNIQ HEALPix pixels
        """
        r = IntervalSet()
        rtmp = IntervalSet()
        last_order = 0
        
        intervals = uniq_is.intervals
        max_res_order = MOC.HPY_MAX_NORDER
        diff_order = max_res_order
        shift_order = 2 * diff_order
        for interval in intervals:
            for j in xrange(interval[0], interval[1]):
                order, ipix = utils.uniq2orderipix(j)
                
                if order != last_order:
                    r = r.union(rtmp)
                    rtmp.clear()
                    last_order = order
                    diff_order = max_res_order - order
                    shift_order = 2 * diff_order

                rtmp.add((ipix << shift_order, (ipix+1) << shift_order))

        r = r.union(rtmp)
        return cls.from_interval_set(r)
    
    @classmethod
    def from_interval_set(cls, interval_set):
        """
        Create a MOC from an IntervalSet (all intervals at deepest norder)
        """
        moc = MOC()
        moc._interval_set = interval_set
        
        return moc
    
    @classmethod
    def from_table(cls, table, ra_column, dec_column, moc_order):
        """
        Create a MOC from a astropy.table.Table
        The user has to specify the columns holding ra and dec (in ICRS)
        """
        if not isinstance(ra_column, str) or not isinstance(dec_column, str):
            raise TypeError

        moc = MOC()
        moc._order = moc_order

        moc.add_position(table[ra_column], table[dec_column], moc_order)
        
        moc._ensure_consistency()

        return moc 
    
    '''
    @classmethod
    def from_coo_list(cls, skycoord_list, max_norder):
        """
        Create a MOC from a list of SkyCoord
        """
        moc = MOC()
        moc._order = max_norder
        
        # very very slow :(, don't really know why
        # Using
        for skycoord in skycoord_list:
            moc.add_position(skycoord.icrs.ra.deg, skycoord.icrs.dec.deg, max_norder)

        return moc
    '''

    @classmethod
    def from_coo_list(cls, skycoords, max_norder):
        """
        Create a MOC from a list of SkyCoord
        """
        if not isinstance(skycoords, SkyCoord):
            raise TypeError

        moc = MOC()
        moc._order = max_norder

        moc.add_position(skycoords.icrs.ra.deg, skycoords.icrs.dec.deg, max_norder)
        moc._ensure_consistency()

        return moc

    def uniq_pixels_iterator(self):
        for uniq_iv in self.to_uniq_interval_set().intervals:
            for uniq in xrange(uniq_iv[0], uniq_iv[1]):
                yield uniq
                
    def best_res_pixels_iterator(self):
        factor = 4**(MOC.HPY_MAX_NORDER - self.max_order)
        for iv in self._interval_set.intervals:
            for val in xrange(iv[0] // factor, iv[1] // factor):
                yield val
    
    def filter_table(self, table, ra_column, dec_column, keep_inside=True):
        """
        Filter an astropy.table.Table to keep only rows inside (or outside) the MOC instance
        Return the (newly created) filtered Table
        """
        filtered_table = Table()
        
        kept_rows = []
        pixels_best_res = set()
        for val in self.best_res_pixels_iterator():
            pixels_best_res.add(val)

        max_order = self.max_order
        nside = 2**max_order
        for row in table:
            theta, phi = utils.radec2thetaphi(row[ra_column], row[dec_column])
            ipix = ang2pix(nside, theta, phi, nest=True)
            if (ipix in pixels_best_res) == keep_inside:
                kept_rows.append(row)

        if len(kept_rows) == 0:
            return Table(names=table.colnames)
        else:
            return Table(rows=kept_rows, names=table.colnames)

    def write(self, path, format='fits', optional_kw_dict = None):
        """
        Serialize a moc in FITS in a given path
        Format can be 'fits' or 'json', though only the fits format is 
        officially supported by the IVOA
        """
        formats = ('fits', 'json')
        if format not in formats:
            raise ValueError('format should be one of %s' % (str(formats)))
        
        self._ensure_consistency()
        uniq_array = []
        for uniq in self.uniq_pixels_iterator():
            uniq_array.append(uniq)
        
        if format == 'fits':
            if self._order is not None:
                moc_order = self._order
            else:
                moc_order = self.max_order
            if moc_order <= 13:
                format = '1J'
            else:
                format = '1K'
            
            tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs([fits.Column(name='UNIQ', format=format, array=np.array(uniq_array))]))
            tbhdu.header['PIXTYPE'] = 'HEALPIX'
            tbhdu.header['ORDERING'] = 'NUNIQ'
            tbhdu.header['COORDSYS'] = 'C'
            tbhdu.header['MOCORDER'] = moc_order
            tbhdu.header['MOCTOOL'] = 'MOCPy'
            if optional_kw_dict:
                for key in optional_kw_dict:
                    tbhdu.header[key] = optional_kw_dict[key]
                    
            thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
            thdulist.writeto(path, clobber=True)
            
        elif format == 'json':
            import json
            
            json_moc = {}
            pix_list = []
            o_order = -1
            for uniq in uniq_array:
                order, ipix = utils.uniq2orderipix(uniq)
                if order!=o_order and o_order>0:
                    json_moc[str(o_order)] = pix_list
                    pix_list = []
                
                o_order = order
                pix_list.append(ipix)
                
            json_moc[str(o_order)] = pix_list
            with open(path, 'w') as h:
                h.write(json.dumps(json_moc, sort_keys=True, indent=2))

        
class MOC_io:
    @staticmethod
    def read_local(path, moc_order=None):
        """
        Read a MOC on the local file system
        """
        return MOC_io.__parse(path, moc_order=moc_order)

    @staticmethod
    def __parse(path_l, moc_order=None):
        result = MOC()

        if not isinstance(path_l, list):
            path_l = [path_l]

        for path in path_l:
            moc = None
            with fits.open(path) as hdulist:
                if isinstance(hdulist[1], fits.hdu.table.BinTableHDU):
                    interval_set = IntervalSet()
                    data = hdulist[1].data.view(np.recarray) # accessing directly recarray dramatically speed up the reading
                    for x in xrange(0, len(hdulist[1].data)):
                        interval_set.add(data[x][0])

                    moc = MOC.from_uniq_interval_set(interval_set)
                elif isinstance(hdulist[1], fits.hdu.ImageHDU):
                    if not moc_order:
                        print('An order must be specified when building a moc from a'
                              'fits image.')
                        raise ValueError
                    # load the image data
                    data_image = fits.getdata(path, ext=0)
                    header = fits.getheader(path)
                    height = header['NAXIS2']
                    width = header['NAXIS1']
                    print(height, width)

                    # use wcs from astropy to locate the image in the world coordinates
                    w = wcs.WCS(header)

                    d_pix = 1
                    X, Y = np.mgrid[0:width:d_pix, 0:height:d_pix]

                    pix_crd = np.dstack((X.ravel(), Y.ravel()))[0]

                    world_pix_crd = w.wcs_pix2world(pix_crd, 1)
                    hp = HEALPix(nside=(1 << moc_order), order='nested', frame=ICRS())

                    print(world_pix_crd)

                    ipix_l = hp.lonlat_to_healpix(world_pix_crd[:, 0] * u.deg,
                                                  world_pix_crd[:, 1] * u.deg)
                    # remove doubles
                    ipix_l = list(set(ipix_l))

                    moc = MOC.from_json({str(moc_order): ipix_l})
                    # the moc computed from the image is not consistent because it is defined
                    # only at the moc order
                    moc._ensure_consistency()

                result = result.union(moc)

        return result

