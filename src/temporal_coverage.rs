//! The `temporal_coverage` module contains methods
//! related to the creation and manipulation of
//! temporal coverages.

use ndarray::{Array, Array1, Array2, Axis};

use crate::coverage;
use pyo3::exceptions;
use pyo3::prelude::PyResult;
use intervals::mocranges::TimeRanges;
use intervals::deser::fits::error::FitsError;

/// Create a temporal coverage from a list of time ranges expressed in jd.
///
/// # Arguments
///
/// * ``min_times`` - The list of inf bounds of the time ranges expressed in **jd**
/// * ``max_times`` - The list of sup bounds of the time ranges expressed in **jd**
///
/// # Errors
///
/// * If the number of ``min_times`` and ``max_times`` do not match.
pub fn from_time_ranges(min_times: Array1<f64>, max_times: Array1<f64>) -> PyResult<Array2<u64>> {
    if min_times.shape() != max_times.shape() {
        Err(exceptions::PyValueError::new_err(
            "min and max ranges have not the same shape",
        ))
    } else {
        if min_times.is_empty() {
            return Ok(Array::zeros((1, 0)));
        }
        let shape = (min_times.shape()[0], 1);

        let min_times = min_times.into_shape(shape).unwrap().mapv(|e| (e * 86400000000_f64) as u64);
        let max_times = max_times.into_shape(shape).unwrap().mapv(|e| (e * 86400000000_f64) as u64);

        let ranges = concatenate![Axis(1), min_times, max_times].to_owned();
        let ranges = coverage::build_time_ranges_from_py(ranges);

        let result: Array2<u64> = ranges.into();
        Ok(result)
    }
}

/// Create a temporal coverage from a list of time ranges expressed in microseconds since jd=0.
///
/// # Arguments
///
/// * ``min_times`` - The list of inf bounds (inclusive) of the time ranges expressed in **microseconds** since **jd=0**
/// * ``max_times`` - The list of sup bounds (exclusive) of the time ranges expressed in **microseconds** since **jd=0**
///
/// # Errors
///
/// * If the number of ``min_times`` and ``max_times`` do not match.
pub fn from_time_ranges_in_microsec_since_jd_origin(min_times: Array1<u64>, max_times: Array1<u64>) -> PyResult<Array2<u64>> {
    if min_times.shape() != max_times.shape() {
        Err(exceptions::PyValueError::new_err(
            "min and max ranges have not the same shape",
        ))
    } else {
        if min_times.is_empty() {
            return Ok(Array::zeros((1, 0)));
        }
        let shape = (min_times.shape()[0], 1);

        let min_times = min_times.into_shape(shape).unwrap();
        let max_times = max_times.into_shape(shape).unwrap();

        let ranges = concatenate![Axis(1), min_times, max_times].to_owned();
        let ranges = coverage::build_time_ranges_from_py(ranges);

        let result: Array2<u64> = ranges.into();
        Ok(result)
    }
}

pub fn to_ascii_str(depth_max: u8, ranges: TimeRanges<u64>) -> String {
    coverage::to_ascii_str(depth_max, ranges)
}
pub fn to_ascii_file(depth_max: u8, ranges: TimeRanges<u64>, path: String) -> std::io::Result<()> {
    coverage::to_ascii_file(depth_max, ranges, path)
}
pub fn to_json_str(depth_max: u8, ranges: TimeRanges<u64>) -> String {
    coverage::to_json_str(depth_max, ranges)
}
pub fn to_json_file(depth_max: u8, ranges: TimeRanges<u64>, path: String) -> std::io::Result<()> {
    coverage::to_json_file(depth_max, ranges, path)
}
pub fn to_fits_file(depth_max: u8, ranges: TimeRanges<u64>, path: String) -> Result<(), FitsError> {
    coverage::to_fits_file(depth_max, ranges, path)
}

pub fn from_ascii_str(ascii: String) -> PyResult<TimeRanges<u64>> {
    coverage::from_ascii_str(ascii).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_ascii_file(path: String) -> PyResult<TimeRanges<u64>> {
    coverage::from_ascii_file(path).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_json_str(json: String) -> PyResult<TimeRanges<u64>> {
    coverage::from_json_str(json).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_json_file(path: String) -> PyResult<TimeRanges<u64>> {
    coverage::from_json_file(path).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_fits_file(path: String) -> PyResult<TimeRanges<u64>> {
    coverage::from_fits_file_time(path).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}