//! The `time_space_coverage` module contains methods
//! related to the creation and manipulation of
//! Time-Space 2D coverages.

use std::convert::TryFrom;

use rayon::prelude::*;

use pyo3::exceptions;
use pyo3::prelude::PyResult;

use intervals::mocqty::{MocQty, Hpx, Time};
use intervals::hpxranges2d::TimeSpaceMoc;
use intervals::mocranges::{HpxRanges, TimeRanges};

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
pub fn create_from_times_positions(
    times: Vec<f64>,
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

        let times = times
            .into_par_iter()
            .map(|t| (t * 86400000000_f64).floor() as u64)
            .collect::<Vec<_>>();

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
pub fn create_from_time_ranges_positions(
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
use pyo3::types::PyList;
use numpy::PyReadonlyArray2;
use pyo3::{ToPyObject, Python};

pub fn from_time_ranges_spatial_coverages(
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

        const ERR_CAST: &'static str = "Cannot cast spatial coverages to Array2<u64>";
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
                .to_owned()
                .into();

            spatial_coverages_res.push(spatial_cov);
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
pub fn contains(
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
        .par_apply(|r, &t, &l, &b| {
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

use ndarray::Array1;
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
    coverage.into()
}

/// Deserialize a Time-Space coverage from FITS
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
pub fn from_fits(data: Array1<i64>) -> PyResult<TimeSpaceMoc<u64, u64>> {
    TimeSpaceMoc::<u64, u64>::try_from(data).map_err(|msg| exceptions::PyValueError::new_err(msg))
}

/// Create a new empty Time-Space coverage
///
/// This method is called in the constructor of the
/// `mocpy.STMOC` class
pub fn new() -> TimeSpaceMoc<u64, u64> {
    // Create new empty coverage
    TimeSpaceMoc::new()
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
pub fn t_min(coverage: &TimeSpaceMoc<u64, u64>) -> PyResult<f64> {
    let t_min = coverage
        .t_min()
        .map_err(|msg| exceptions::PyValueError::new_err(msg))?;

    Ok((t_min as f64) / 86400000000_f64)
}

/// Returns the maximum value along the `T` dimension
///
/// # Errors
///
/// When the `TimeSpaceMoc<T, S>` is empty.
pub fn t_max(coverage: &TimeSpaceMoc<u64, u64>) -> PyResult<f64> {
    let t_max = coverage
        .t_max()
        .map_err(|msg| exceptions::PyValueError::new_err(msg))?;

    Ok((t_max as f64) / 86400000000_f64)
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
