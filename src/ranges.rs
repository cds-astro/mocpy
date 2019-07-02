use intervals::NestedRanges;
use intervals::bounded::Bounded;
use std::ops::Range;

use rayon::prelude::*;

pub fn create_from_position(lon: Vec<f64>, lat: Vec<f64>, depth: u8) -> NestedRanges<u64> {
    let mut data = Vec::<Range<u64>>::with_capacity(lon.len());
    data.resize(lon.len(), 0..1);

    let shift = (<u64>::MAXDEPTH as u8 - depth) << 1;
    let layer = healpix::nested::get_or_create(depth);

    data.par_iter_mut()
        .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
        .for_each(|(p, (l, b))| {
            let pix = layer.hash(l, b);
            
            let e1 = pix << shift;
            let e2 = (pix + 1) << shift;
            *p = e1..e2;
        });

    NestedRanges::<u64>::new(data, None, true)
}