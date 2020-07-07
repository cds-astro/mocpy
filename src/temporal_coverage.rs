//! The `temporal_coverage` module contains methods
//! related to the creation and manipulation of
//! temporal coverages.

use ndarray::{Array, Array1, Array2, Axis};

use crate::coverage;
use pyo3::exceptions;
use pyo3::prelude::PyResult;
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
        Err(exceptions::ValueError::py_err(
            "min and max ranges have not the same shape",
        ))
    } else {
        if min_times.is_empty() {
            return Ok(Array::zeros((1, 0)));
        }
        let shape = (min_times.shape()[0], 1);

        let min_times = min_times.into_shape(shape).unwrap();
        let max_times = max_times.into_shape(shape).unwrap();

        let min_times = &min_times * &Array::from_elem(shape, 86400000000_f64);
        let min_times = min_times.mapv(|e| e as u64);

        let max_times = &max_times * &Array::from_elem(shape, 86400000000_f64);
        let max_times = max_times.mapv(|e| e as u64);

        let ranges = stack![Axis(1), min_times, max_times].to_owned();
        let ranges = coverage::create_nested_ranges_from_py(ranges).make_consistent();

        let result: Array2<u64> = ranges.into();
        Ok(result)
    }
}
