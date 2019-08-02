use intervals::nestedranges::NestedRanges;

use intervals::bounded::Bounded;
use std::ops::Range;

use rayon::prelude::*;
use pyo3::exceptions;
use pyo3::prelude::PyResult;

/// Create a spatial coverage from a list of sky coordinates
/// 
/// # Arguments
/// 
/// * ``lon`` - The longitudes of the sky coordinates.
/// * ``lat`` - The latitudes of the sky coordinates.
/// * ``depth`` - The depth at which HEALPix cell indices
///   will be computed. This will correspond to the depth
///   of the coverage once it will be created.
/// 
/// # Precondition
/// 
/// ``lon`` and ``lat`` are expressed in radians.
/// They are valid because they come from
/// `astropy.coordinates.Quantity` objects.
/// 
/// # Errors
/// 
/// * If the number of longitudes and latitudes do not match.
/// * If ``depth`` is not in `[0, <u64>::MAXDEPTH] = [0, 29]`
pub fn create_from_position(lon: Vec<f64>, lat: Vec<f64>, depth: i8) -> PyResult<NestedRanges<u64>> {
    if lon.len() != lat.len() {
        return Err(exceptions::ValueError::py_err("Longitudes and Latitudes \
            do not have the same shapes."));
    }
    if depth < 0 || depth > <u64>::MAXDEPTH {
        return Err(exceptions::ValueError::py_err(
            format!("Depth must be comprised between in [0, {0}]", <u64>::MAXDEPTH)
        ));
    }

    let mut data = Vec::<Range<u64>>::with_capacity(lon.len());
    data.resize(lon.len(), 0..1);

    let shift = (<u64>::MAXDEPTH - depth) << 1;
    let layer = healpix::nested::get_or_create(depth as u8);

    data.par_iter_mut()
        .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
        .for_each(|(p, (l, b))| {
            let pix = layer.hash(l, b);
            
            let e1 = pix << shift;
            let e2 = (pix + 1) << shift;
            *p = e1..e2;
        });

    ;
    let result = NestedRanges::<u64>::new(data).make_consistent();
    Ok(result)
}

use intervals::nestedranges2d::NestedRanges2D;
/// Create a time-spatial coverage (2D) from a list of sky coordinates
/// and times.
/// 
/// # Arguments
/// 
/// * ``times`` - The times expressed in jd.
/// * ``lon`` - The longitudes of the sky coordinates.
/// * ``lat`` - The latitudes of the sky coordinates.
/// * ``dt`` - The depth along the time (i.e. `T`) axis.
/// * ``ds`` - The depth at which HEALPix cell indices
///   will be computed.
/// 
/// # Precondition
/// 
/// * ``lon`` and ``lat`` are expressed in radians.
/// They are valid because they come from
/// `astropy.coordinates.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
/// 
/// # Errors
/// 
/// If the number of longitudes, latitudes and times do not match.
pub fn create_from_time_position(times: Vec<f64>, lon: Vec<f64>, lat: Vec<f64>, dt: i8, ds: i8) -> PyResult<NestedRanges2D<u64, u64>> {
    if dt < 0 || dt > u64::MAXDEPTH {
        Err(exceptions::ValueError::py_err(
            format!("Time depth must be in [0, {0}]", <u64>::MAXDEPTH)
        ))
    } else if ds < 0 || ds > u64::MAXDEPTH {
        Err(exceptions::ValueError::py_err(
            format!("Space depth must be in [0, {0}]", <u64>::MAXDEPTH)
        ))
    } else {
        if times.len() != lon.len() || times.len() != lat.len() {
            return Err(exceptions::ValueError::py_err(
                "Times, longitudes and latitudes do not have
                 the same shapes."));
        }

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

/// Degrade a spatial coverage to a specific depth.
/// 
/// # Arguments
/// 
/// * ``ranges`` - The spatial coverage ranges to degrade.
/// * ``depth`` - The depth to degrade the spatial coverage to.
/// 
/// # Errors
/// 
/// * ``depth`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 29]`
pub fn degrade_nested_ranges(ranges: &mut NestedRanges<u64>, depth: i8) -> PyResult<()> {
    if depth < 0 || depth > <u64>::MAXDEPTH {
        Err(exceptions::ValueError::py_err(
            format!("Depth must be in [0, {0}]", <u64>::MAXDEPTH)
        ))
    } else {
        ranges.degrade(depth);
        Ok(())
    }
}

use ndarray::Array2;
/// Cast an `Array2<u64>` coming from MOCPy python code to
/// a `NestedRanges<u64>` object.
pub fn create_nested_ranges_from_py(data: Array2<u64>) -> NestedRanges<u64> {
    data.into()
}

use intervals::uniqranges::UniqRanges;
/// Cast an `Array2<u64>` coming from MOCPy python code to
/// an `UniqRanges<u64>` object.
pub fn create_uniq_ranges_from_py(data: Array2<u64>) -> UniqRanges<u64> {
    data.into()
}