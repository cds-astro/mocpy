//! The `time_space_coverage` module contains methods
//! related to the creation and manipulation of
//! Time-Space 2D coverages.

use std::path::Path;
use std::fs;
use std::fs::{File};
use std::io::{BufReader, BufWriter};
use std::error::Error;

use rayon::prelude::*;

use ndarray::Array1;

use numpy::PyReadonlyArray2;

use pyo3::types::PyList;
use pyo3::{ToPyObject, Python};
use pyo3::exceptions;
use pyo3::prelude::PyResult;


use moc::deser::fits::{from_fits_ivoa, MocIdxType, MocQtyType, ranges2d_to_fits_ivoa, STMocType};
use moc::deser::ascii::{AsciiError, moc2d_from_ascii_ivoa};
use moc::deser::json::cellmoc2d_from_json_aladin;
use moc::qty::{MocQty, Hpx, Time};
use moc::hpxranges2d::TimeSpaceMoc;
use moc::elemset::range::{HpxRanges, TimeRanges};
use moc::moc2d::{
    CellMOC2Iterator, CellMOC2IntoIterator,
    CellOrCellRangeMOC2Iterator, CellOrCellRangeMOC2IntoIterator
};

use crate::ndarray_fromto::{
    array2_to_mocranges, 
    hpxranges2d_to_array1_i64, 
    try_array1_i64_to_hpx_ranges2, try_array1_u64_to_hpx_ranges2
};

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
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
/// 
/// # Remark 
/// 
/// Method kept temporarily to ensure backward compatibility.
/// 
pub fn create_from_times_positions_approx(
    times: Vec<f64>,
    lon: Vec<f64>,
    lat: Vec<f64>,
    dt: u8,
    ds: u8,
) -> PyResult<TimeSpaceMoc<u64, u64>> {
    let times = times
      .into_par_iter()
      .map(|t| (t * 86400000000_f64).floor() as u64)
      .collect::<Vec<_>>();
    create_from_times_positions(times, lon, lat, dt, ds)
}

/// Create a time-spatial coverage (2D) from a list of sky coordinates
/// and times.
///
/// # Arguments
///
/// * ``times`` - The times expressed in microsecond since jd=0.
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
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
pub fn create_from_times_positions(
    times: Vec<u64>,
    lon: Vec<f64>,
    lat: Vec<f64>,
    dt: u8,
    ds: u8,
) -> PyResult<TimeSpaceMoc<u64, u64>> {
    if dt > Time::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Time depth must be in [0, {0}]",
            Time::<u64>::MAX_DEPTH
        )))
    } else if ds > Hpx::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Space depth must be in [0, {0}]",
            Hpx::<u64>::MAX_DEPTH
        )))
    } else {
        if times.len() != lon.len() || times.len() != lat.len() {
            return Err(exceptions::PyValueError::new_err(
                "Times, longitudes and latitudes do not have
                 the same shapes.",
            ));
        }

        let mut ipix = vec![0; lon.len()];

        let layer = healpix::nested::get(ds as u8);
        ipix.par_iter_mut()
          .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
          .for_each(|(p, (l, b))| {
              *p = layer.hash(l, b);
          });

        Ok(TimeSpaceMoc::<u64, u64>::create_from_times_positions(
            times, ipix, dt, ds,
        ))
    }
}

/// Create a time-spatial coverage (2D) from a list of sky coordinates
/// and ranges of times.
///
/// # Arguments
///
/// * ``times_start`` - The starting times expressed in jd.
/// * ``times_end`` - The ending times expressed in jd.
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
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
/// 
/// # Remark 
/// 
/// Method kept temporarily to ensure backward compatibility.
pub fn create_from_time_ranges_positions_approx(
    times_start: Vec<f64>,
    times_end: Vec<f64>,
    dt: u8,
    lon: Vec<f64>,
    lat: Vec<f64>,
    ds: u8,
) -> PyResult<TimeSpaceMoc<u64, u64>> {
    if ds > Hpx::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Space depth must be in [0, {0}]",
            Hpx::<u64>::MAX_DEPTH
        )))
    } else if dt > Time::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Time depth must be in [0, {0}]",
            Time::<u64>::MAX_DEPTH
        )))
    } else {
        if times_start.len() != lon.len() ||
           times_start.len() != lat.len() ||
           times_start.len() != times_end.len() {
            return Err(exceptions::PyValueError::new_err(
                "Times, longitudes and latitudes do not have the same shapes.",
            ));
        }

        let mut ipix = vec![0; lon.len()];

        let layer = healpix::nested::get(ds as u8);
        ipix.par_iter_mut()
            .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
            .for_each(|(p, (l, b))| {
                *p = layer.hash(l, b);
            });

        let times = times_start.into_par_iter()
            .zip_eq(times_end.into_par_iter())
            .filter_map(|(t1, t2)| {
                let t1 = (t1 * 86400000000_f64).floor() as u64;
                let t2 = (t2 * 86400000000_f64).floor() as u64;
                if t2 > t1 {
                    Some(t1..t2)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        if times.len() != ipix.len() {
            return Err(exceptions::PyValueError::new_err(
                "Number of time ranges and sky coordinates do not match.",
            ));
        }

        Ok(TimeSpaceMoc::<u64, u64>::create_from_time_ranges_positions(
            times, ipix, dt, ds,
        ))
    }
}

/// Create a time-spatial coverage (2D) from a list of sky coordinates
/// and ranges of times.
///
/// # Arguments
///
/// * ``times_start`` - The starting times expressed in microseconds since jd=0.
/// * ``times_end`` - The ending times expressed in microseconds since jd=0.
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
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
pub fn create_from_time_ranges_positions(
    times_start: Vec<u64>,
    times_end: Vec<u64>,
    dt: u8,
    lon: Vec<f64>,
    lat: Vec<f64>,
    ds: u8,
) -> PyResult<TimeSpaceMoc<u64, u64>> {
    if ds > Hpx::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Space depth must be in [0, {0}]",
            Hpx::<u64>::MAX_DEPTH
        )))
    } else if dt > Time::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Time depth must be in [0, {0}]",
            Time::<u64>::MAX_DEPTH
        )))
    } else {
        if times_start.len() != lon.len() ||
          times_start.len() != lat.len() ||
          times_start.len() != times_end.len() {
            return Err(exceptions::PyValueError::new_err(
                "Times, longitudes and latitudes do not have the same shapes.",
            ));
        }

        let mut ipix = vec![0; lon.len()];

        let layer = healpix::nested::get(ds as u8);
        ipix.par_iter_mut()
          .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
          .for_each(|(p, (l, b))| {
              *p = layer.hash(l, b);
          });

        let times = times_start.into_par_iter()
          .zip_eq(times_end.into_par_iter())
          .map(|(t1, t2)| t1..t2)
          .collect::<Vec<_>>();

        if times.len() != ipix.len() {
            return Err(exceptions::PyValueError::new_err(
                "Number of time ranges and sky coordinates do not match.",
            ));
        }

        Ok(TimeSpaceMoc::<u64, u64>::create_from_time_ranges_positions(
            times, ipix, dt, ds,
        ))
    }
}




/// Create a time-spatial coverage (2D) from a list of cones
/// and time ranges.
///
/// # Arguments
///
/// * ``times_start`` - The starting times expressed in jd.
/// * ``times_end`` - The ending times expressed in jd.
/// * ``lon`` - The longitudes of the sky coordinates.
/// * ``lat`` - The latitudes of the sky coordinates.
/// * ``radius`` - The radiuses of the cones.
/// * ``dt`` - The depth along the time (i.e. `T`) axis.
/// * ``ds`` - The depth at which HEALPix cell indices
///   will be computed.
///
/// # Precondition
///
/// * ``lon`` and ``lat`` are expressed in radians.
/// They are valid because they come from
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
/// 
/// # Remark 
/// 
/// Method kept temporarily to ensure backward compatibility.
/// 
pub fn from_time_ranges_spatial_coverages_approx(
    py: Python,
    times_start: Vec<f64>,
    times_end: Vec<f64>,
    dt: u8,
    spatial_coverages: &PyList,
) -> PyResult<TimeSpaceMoc<u64, u64>> {
    if dt > Time::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Time depth must be in [0, {0}]",
            Time::<u64>::MAX_DEPTH
        )))
    } else {
        if times_start.len() != times_end.len() {
            return Err(exceptions::PyValueError::new_err(
                "Invalid times.",
            ));
        }

        const ERR_CAST: &str = "Cannot cast spatial coverages to Array2<u64>";
        let mut spatial_coverages_res: Vec<HpxRanges<u64>> = vec![];

        for spatial_cov in spatial_coverages.into_iter() {
            // Loop over the python list and cast its elements
            // to numpy PyArray2<u64> types first, to ndarray Array2<u64> secondly
            // and finally to HpxRanges<u64>
            let spatial_cov = spatial_cov
                .to_object(py)
                .extract::<PyReadonlyArray2<u64>>(py)
                .map_err(|_| exceptions::PyValueError::new_err(ERR_CAST))?
                .as_array()
                .to_owned();

            spatial_coverages_res.push(array2_to_mocranges(spatial_cov));
        }

        let times = times_start.into_par_iter()
            .zip_eq(times_end.into_par_iter())
            .filter_map(|(t1, t2)| {
                let t1 = (t1 * 86400000000_f64).floor() as u64;
                let t2 = (t2 * 86400000000_f64).floor() as u64;
                if t2 > t1 {
                    Some(t1..t2)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        Ok(TimeSpaceMoc::<u64, u64>::create_from_time_ranges_spatial_coverage(
            times, spatial_coverages_res, dt,
        ))
    }
}

/// Create a time-spatial coverage (2D) from a list of cones
/// and time ranges.
///
/// # Arguments
///
/// * ``times_start`` - The starting times expressed in microseconds since jd=0.
/// * ``times_end`` - The ending times expressed in  microseconds since jd=0.
/// * ``lon`` - The longitudes of the sky coordinates.
/// * ``lat`` - The latitudes of the sky coordinates.
/// * ``radius`` - The radiuses of the cones.
/// * ``dt`` - The depth along the time (i.e. `T`) axis.
/// * ``ds`` - The depth at which HEALPix cell indices
///   will be computed.
///
/// # Precondition
///
/// * ``lon`` and ``lat`` are expressed in radians.
/// They are valid because they come from
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
pub fn from_time_ranges_spatial_coverages(
    py: Python,
    times_start: Vec<u64>,
    times_end: Vec<u64>,
    dt: u8,
    spatial_coverages: &PyList,
) -> PyResult<TimeSpaceMoc<u64, u64>> {
    if dt > Time::<u64>::MAX_DEPTH {
        Err(exceptions::PyValueError::new_err(format!(
            "Time depth must be in [0, {0}]",
            Time::<u64>::MAX_DEPTH
        )))
    } else {
        if times_start.len() != times_end.len() {
            return Err(exceptions::PyValueError::new_err(
                "Invalid times.",
            ));
        }

        const ERR_CAST: &str = "Cannot cast spatial coverages to Array2<u64>";
        let mut spatial_coverages_res: Vec<HpxRanges<u64>> = vec![];

        for spatial_cov in spatial_coverages.into_iter() {
            // Loop over the python list and cast its elements
            // to numpy PyArray2<u64> types first, to ndarray Array2<u64> secondly
            // and finally to HpxRanges<u64>
            let spatial_cov = spatial_cov
              .to_object(py)
              .extract::<PyReadonlyArray2<u64>>(py)
              .map_err(|_| exceptions::PyValueError::new_err(ERR_CAST))?
              .as_array()
              .to_owned();

            spatial_coverages_res.push(array2_to_mocranges(spatial_cov));
        }

        let times = times_start.into_par_iter()
          .zip_eq(times_end.into_par_iter())
          .map(|(t1, t2)| t1..t2)
          .collect::<Vec<_>>();

        Ok(TimeSpaceMoc::<u64, u64>::create_from_time_ranges_spatial_coverage(
            times, spatial_coverages_res, dt,
        ))
    }
}

use ndarray::Zip;
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
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
/// # Remark 
/// 
/// Method kept temporarily to ensure backward compatibility.
/// 
pub fn contains_approx(
    coverage: &TimeSpaceMoc<u64, u64>,
    time: Array1<f64>,
    lon: Array1<f64>,
    lat: Array1<f64>,
    result: &mut Array1<bool>,
) -> PyResult<()> {
    if time.len() != lon.len() || time.len() != lat.len() {
        return Err(exceptions::PyValueError::new_err(
            "Times, longitudes and latitudes do not have the same shapes.",
        ));
    }

    // Retrieve the spatial depth of the Time-Space coverage
    let (_, s_depth) = coverage.compute_min_depth();
    let layer = healpix::nested::get(s_depth as u8);
    let shift = (Hpx::<u64>::MAX_DEPTH - s_depth) << 1;
    Zip::from(result)
        .and(&time)
        .and(&lon)
        .and(&lat)
        .par_for_each(|r, &t, &l, &b| {
            // Compute the HEALPix cell range at the max depth
            // along the spatial dimension
            let pix = layer.hash(l, b);
            let e1 = pix << shift;
            let e2 = (pix + 1) << shift;

            // Convert the observation time in µs
            let t = (t * 86400000000_f64).floor() as u64;
            // Check whether the (time in µs, HEALPix cell nested range)
            // is contained into the Spatial-Time coverage
            *r = coverage.contains(t, &(e1..e2));
        });

    Ok(())
}

/// Create a time-spatial coverage (2D) from a list of sky coordinates
/// and times.
///
/// # Arguments
///
/// * ``times`` - The times expressed in microseconds since jd=0.
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
/// `astropy.units.Quantity` objects.
/// * ``times`` are expressed in jd and are coming
/// from `astropy.time.Time` objects.
///
/// # Errors
///
/// If the number of longitudes, latitudes and times do not match.
pub fn contains(
    coverage: &TimeSpaceMoc<u64, u64>,
    time: Array1<u64>,
    lon: Array1<f64>,
    lat: Array1<f64>,
    result: &mut Array1<bool>,
) -> PyResult<()> {
    if time.len() != lon.len() || time.len() != lat.len() {
        return Err(exceptions::PyValueError::new_err(
            "Times, longitudes and latitudes do not have the same shapes.",
        ));
    }

    // Retrieve the spatial depth of the Time-Space coverage
    let (_, s_depth) = coverage.compute_min_depth();
    let layer = healpix::nested::get(s_depth as u8);
    let shift = (Hpx::<u64>::MAX_DEPTH - s_depth) << 1;
    Zip::from(result)
      .and(&time)
      .and(&lon)
      .and(&lat)
      .par_for_each(|r, &t, &l, &b| {
          // Compute the HEALPix cell range at the max depth
          // along the spatial dimension
          let pix = layer.hash(l, b);
          let e1 = pix << shift;
          let e2 = (pix + 1) << shift;
          
          // Check whether the (time in µs, HEALPix cell nested range)
          // is contained into the Spatial-Time coverage
          *r = coverage.contains(t, &(e1..e2));
      });

    Ok(())
}


/// Returns the union of the ranges along the `T` axis for which their
/// `S` ranges is contained in ``y``
///
/// # Arguments
///
/// * ``y``- The set of ranges along the `S` axis.
/// * ``coverage`` - The input coverage.
///
/// # Algorithm
///
/// This method checks for all the `S` axis ranges of ``coverage`` that
/// lie into the range set ``y``.
///
/// It then performs the union of the `T` axis ranges corresponding to the
/// matching ranges along the `S` axis.
pub fn project_on_first_dim(
    y: &HpxRanges<u64>,
    coverage: &TimeSpaceMoc<u64, u64>,
) -> TimeRanges<u64> {
    TimeSpaceMoc::project_on_first_dim(y, coverage)
}

/// Returns the union of the ranges along the `S` axis for which their
/// `T` ranges intersect ``x``
///
/// # Arguments
///
/// * ``x``- The set of ranges along the `T` axis.
/// * ``coverage`` - The input coverage
///
/// # Algorithm
///
/// This method checks for all the `T` axis ranges of ``coverage`` that
/// lie into the range set ``x``.
///
/// It then performs the union of the `S` axis ranges corresponding to the
/// matching ranges along the `T` axis.
pub fn project_on_second_dim(
    x: &TimeRanges<u64>,
    coverage: &TimeSpaceMoc<u64, u64>,
) -> HpxRanges<u64> {
    TimeSpaceMoc::project_on_second_dim(x, coverage)
}

/// Create a Array1<i64> from a TimeSpaceMoc<u64, u64>
///
/// This is used when storing a STMOC into a FITS file
///
/// # Info
///
/// The output Array1 stores the STMOC under the nested format.
/// Its memory layout contains each time range followed by the
/// list of space ranges referred to that time range.
/// Time ranges are negatives so that one can distinguish them
/// from space ranges.
///
/// Content example of an Array1 coming from a FITS file:
/// int64[] = {-1, -3, 3, 5, 10, 12, 13, 18, -5, -6, 0, 1}
pub fn to_fits(coverage: &TimeSpaceMoc<u64, u64>) -> Array1<i64> {
    hpxranges2d_to_array1_i64(coverage)
}

/// Deserialize a Time-Space coverage from FITS, using the pre-version 2.0 MOC standard.
///
/// # Context
///
/// This is wrapped around the `from_fits` method
/// of MOCPy to load a Time-Space coverage from a
/// FITS file.
///
/// # Arguments
///
/// * ``data`` - A 1d array buffer containing the time and
///   space axis ranges data.
///
/// # Errors
///
/// The `Array1` object stores the Time-Space coverage
/// under the nested format.
/// Its memory layout contains each time range followed by the
/// list of space ranges referred to that time range.
/// Time ranges are negatives so that one can distinguish them
/// from space ranges.
///
/// This method returns a `PyValueError` if the `Array1` is not
/// defined as above.
pub fn from_fits_pre_v2(data: Array1<i64>) -> PyResult<TimeSpaceMoc<u64, u64>> {
    try_array1_i64_to_hpx_ranges2::<Time<u64>>(data).map_err(exceptions::PyValueError::new_err)
}

/// Deserialize a Time-Space coverage from FITS, using the MOC2.0 standard.
///
/// # Context
///
/// This is wrapped around the `from_fits` method
/// of MOCPy to load a Time-Space coverage from a
/// FITS file.
///
/// # Arguments
///
/// * ``data`` - A 1d array buffer containing the time and
///   space axis ranges data.
///
/// # Errors
///
/// The `Array1` object stores the Time-Space coverage
/// under the nested format.
/// Its memory layout contains each time range followed by the
/// list of space ranges referred to that time range.
/// Time ranges are negatives so that one can distinguish them
/// from space ranges.
///
/// This method returns a `PyValueError` if the `Array1` is not
/// defined as above.
pub fn from_fits(data: Array1<u64>) -> PyResult<TimeSpaceMoc<u64, u64>> {
    try_array1_u64_to_hpx_ranges2::<Time<u64>>(data).map_err(exceptions::PyValueError::new_err)
}

/// Deserialize a Time-Space coverage from a FITS file, using the MOC2.0 standard.
///
/// # Arguments
///
/// * ``path`` - path of the ST-MOC fits file
///
/// # Warning
/// 
/// This function is not compatible with pre-v2.0 ST-MOCs.
/// 
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in reading the FITS file
/// (I/O error, format not recognized, ...).
pub fn from_fits_file(path: &Path) -> PyResult<TimeSpaceMoc<u64, u64>> {
    // See https://github.com/PyO3/pyo3/blob/88d86a65aa78bf2d001753f994fe3f6db1d8d75e/src/err/impls.rs
    let file = File::open(&path).map_err(exceptions::PyValueError::new_err)?;
    let reader = BufReader::new(file);
    match from_fits_ivoa(reader).map_err(|err| exceptions::PyIOError::new_err(err.to_string()))? {
        MocIdxType::U64(MocQtyType::TimeHpx(STMocType::V2(it))) => Ok(TimeSpaceMoc::<u64, u64>::from_ranges_it(it)),
        MocIdxType::U64(MocQtyType::TimeHpx(STMocType::PreV2(it))) => Ok(TimeSpaceMoc::<u64, u64>::from_ranges_it(it)),
        _ => Err(exceptions::PyIOError::new_err("Only ST-MOC of u64 ranges supported!")),
    }
}

/// Deserialize a Time-Space coverage from a JSON string, using the MOC2.0 standard.
///
/// # Arguments
///
/// * ``json`` - the ST-MOC JSON string
///
pub fn from_json_str(json: String) -> PyResult<TimeSpaceMoc<u64, u64>> {
    let cellmoc2 = cellmoc2d_from_json_aladin::<u64, Time::<u64>, u64, Hpx::<u64>>(&json)
      .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))?;
    Ok(TimeSpaceMoc::<u64, u64>::from(cellmoc2.into_cell_moc2_iter()))
}

/// Deserialize a Time-Space coverage from a JSON file.
///
/// # Arguments
///
/// * ``path`` - path of the ST-MOC JSON file
///
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in reading the JSON file
/// (I/O error, format not recognized, ...).
pub fn from_json_file(path: &Path) -> PyResult<TimeSpaceMoc<u64, u64>> {
    let file_content = fs::read_to_string(&path)
      .map_err(exceptions::PyIOError::new_err)?;
    from_json_str(file_content)
}

/// Deserialize a Time-Space coverage from a ASCII string, using the MOC2.0 standard.
///
/// # Arguments
///
/// * ``ascii`` - the ST-MOC ASCIi string
///
pub fn from_ascii_str(ascii: String) -> PyResult<TimeSpaceMoc<u64, u64>> {
    let cellmoc2 = moc2d_from_ascii_ivoa::<u64, Time::<u64>, u64, Hpx::<u64>>(&ascii)
      .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))?;
    Ok(TimeSpaceMoc::<u64, u64>::from(cellmoc2.into_cellcellrange_moc2_iter()))
}

/// Deserialize a Time-Space coverage from a ASCII file, using the MOC2.0 standard.
///
/// # Arguments
///
/// * ``path`` - path of the ST-MOC ASCII file
///
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in reading the ASCII file
/// (I/O error, format not recognized, ...).
pub fn from_ascii_file(path: &Path) -> PyResult<TimeSpaceMoc<u64, u64>> {
    let file_content = fs::read_to_string(&path)
      .map_err(exceptions::PyIOError::new_err)?;
    from_ascii_str(file_content)
}


/// Create a new empty Time-Space coverage
///
/// This method is called in the constructor of the
/// `mocpy.STMOC` class
pub fn new() -> TimeSpaceMoc<u64, u64> {
    // Create new empty coverage
    Default::default()
}

/// Compute the depth of the coverage
///
/// # Returns
///
/// A tuple containing two values:
///
/// * The maximum depth along the `T` axis
/// * The maximum depth along the `S` axis
///
/// # Info
///
/// If the `TimeSpaceMoc<T, S>` is empty, the depth returned
/// is set to (0, 0)
pub fn depth(coverage: &TimeSpaceMoc<u64, u64>) -> (u8, u8) {
    coverage.compute_min_depth()
}

/// Returns the minimum value along the `T` dimension
///
/// # Errors
///
/// When the `TimeSpaceMoc<T, S>` is empty.
///
/// # Remark 
/// 
/// Method kept temporarily to ensure backward compatibility.
/// 
pub fn t_min_jd(coverage: &TimeSpaceMoc<u64, u64>) -> PyResult<f64> {
    let t_min = coverage.t_min().map_err(exceptions::PyValueError::new_err)?;

    Ok((t_min as f64) / 86400000000_f64)
}

/// Returns the maximum value along the `T` dimension
///
/// # Errors
///
/// When the `TimeSpaceMoc<T, S>` is empty.
/// 
/// # Remark 
/// 
/// Method kept temporarily to ensure backward compatibility.
/// 
pub fn t_max_jd(coverage: &TimeSpaceMoc<u64, u64>) -> PyResult<f64> {
    let t_max = coverage.t_max().map_err(exceptions::PyValueError::new_err)?;

    Ok((t_max as f64) / 86400000000_f64)
}

/// Returns the minimum value along the `T` dimension
///
/// # Errors
///
/// When the `TimeSpaceMoc<T, S>` is empty.
pub fn t_min_mircosecond_since_jd_org(coverage: &TimeSpaceMoc<u64, u64>) -> PyResult<u64> {
    coverage.t_min().map_err(exceptions::PyValueError::new_err)
}

/// Returns the maximum value along the `T` dimension
///
/// # Errors
///
/// When the `TimeSpaceMoc<T, S>` is empty.
pub fn t_max_mircosecond_since_jd_org(coverage: &TimeSpaceMoc<u64, u64>) -> PyResult<u64> {
    coverage.t_max().map_err(exceptions::PyValueError::new_err)
}

/// Performs the union between two `TimeSpaceMoc<T, S>`
///
/// # Arguments
///
/// * ``other`` - The other `TimeSpaceMoc<T, S>` to
///   perform the union with.
pub fn union(
    coverage_left: &TimeSpaceMoc<u64, u64>,
    coverage_right: &TimeSpaceMoc<u64, u64>,
) -> TimeSpaceMoc<u64, u64> {
    coverage_left.union(coverage_right)
}

/// Performs the intersection between two `TimeSpaceMoc<T, S>`
///
/// # Arguments
///
/// * ``other`` - The other `TimeSpaceMoc<T, S>` to
///   perform the intersection with.
pub fn intersection(
    coverage_left: &TimeSpaceMoc<u64, u64>,
    coverage_right: &TimeSpaceMoc<u64, u64>,
) -> TimeSpaceMoc<u64, u64> {
    coverage_left.intersection(coverage_right)
}

/// Performs the difference between two `TimeSpaceMoc<T, S>`
///
/// # Arguments
///
/// * ``other`` - The other `TimeSpaceMoc<T, S>` to
///   perform the difference with.
pub fn difference(
    coverage_left: &TimeSpaceMoc<u64, u64>,
    coverage_right: &TimeSpaceMoc<u64, u64>,
) -> TimeSpaceMoc<u64, u64> {
    coverage_left.difference(coverage_right)
}


pub fn to_ascii_str(depth_max_t: u8, depth_max_s: u8, coverage: &TimeSpaceMoc<u64, u64>) -> String {
    let mut ascii = Vec::new();
    coverage.time_space_iter(depth_max_t, depth_max_s)
      .into_cellcellrange_moc2_iter()
      .to_ascii_ivoa(Some(80), false, &mut ascii)
      .unwrap(); // unwrap since we do not write in a file but in a string
    unsafe{ String::from_utf8_unchecked(ascii) }
}

pub fn to_ascii_file(depth_max_t: u8, depth_max_s: u8, coverage: &TimeSpaceMoc<u64, u64>, path: String) ->Result<(), AsciiError> {
    let file = File::create(Path::new(&path))?;
    let writer = BufWriter::new(file);
    coverage.time_space_iter(depth_max_t, depth_max_s)
      .into_cellcellrange_moc2_iter()
      .to_ascii_ivoa(Some(80), false, writer)
}

pub fn to_json_str(depth_max_t: u8, depth_max_s: u8, coverage: &TimeSpaceMoc<u64, u64>) -> String {
    let mut json = Vec::new();
    coverage.time_space_iter(depth_max_t, depth_max_s)
      .into_cell_moc2_iter()
      .to_json_aladin(&Some(80), &mut json)
      .unwrap(); // unwrap since we do not write in a file but in a string
    unsafe{ String::from_utf8_unchecked(json) }
}

pub fn to_json_file(depth_max_t: u8, depth_max_s: u8, coverage: &TimeSpaceMoc<u64, u64>, path: String) -> std::io::Result<()> {
    let file = File::create(Path::new(&path))?;
    let writer = BufWriter::new(file);
    coverage.time_space_iter(depth_max_t, depth_max_s)
      .into_cell_moc2_iter()
      .to_json_aladin(&Some(80), writer)
}


/// Write the given coverage into a FITS file of given path.
// TODO: add MOCID and MOCTYPE using Option<&PyDic>
pub fn to_fits_file(depth_max_t: u8, depth_max_s: u8, coverage: &TimeSpaceMoc<u64, u64>, path: &Path) -> Result<(), Box<dyn Error>> {
    let file = File::create(path).map_err(Box::new)?;
    let writer = BufWriter::new(file);
    /*coverage.time_space_iter(depth_max_t, depth_max_s)
      .to_fits_ivoa(None, None, writer)
      .map_err(|err| exceptions::PyIOError::new_err(err.to_string()))*/
    ranges2d_to_fits_ivoa(
        coverage.time_space_iter(depth_max_t, depth_max_s),
        None, None, writer
    ).map_err(|e| e.into())
}
