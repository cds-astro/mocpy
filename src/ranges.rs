use intervals::NestedRanges;
use intervals::ranges2d::NestedRanges2D;
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

pub fn create_from_time_position(times: Vec<f64>, lon: Vec<f64>, lat: Vec<f64>, dt: i8, ds: i8) -> Result<NestedRanges2D<u64, u64>, &'static str> {
    if dt > u64::MAXDEPTH {
        Err("Time depth must be <= 29")
    } else if ds > u64::MAXDEPTH {
        Err("Space depth must be <= 29")
    } else {
        let mut ipix = vec![0; lon.len()];

        let layer = healpix::nested::get_or_create(ds as u8);
        ipix.par_iter_mut()
            .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
            .for_each(|(p, (l, b))| {
                *p = layer.hash(l, b);
            });
        
        let times = times.into_par_iter()
            .map(|t| {
                (t * 86400000000_f64).floor() as u64
            })
            .collect::<Vec<_>>();

        Ok(NestedRanges2D::<u64, u64>::create_quantity_space_coverage(times, ipix, dt, ds))
        //Ok(NestedRanges2D::<u64, u64>::create_space_quantity_coverage(ipix, times, ds as i8, dt as i8))
    }
}

/*
use ndarray::{Array2, Array3};
pub fn convert_nested_ranges_2d(ranges: NestedRanges2D<u64, u64>) -> (Array2<u64>, Array3<u64>) {
    let x = NestedRanges::<u64>::new(ranges.ranges.x, None, false);
}*/