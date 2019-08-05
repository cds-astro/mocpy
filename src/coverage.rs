//! The `coverage` module contains general methods
//! concerning 1d coverages.
//! 
//! These are not specific to a particular dimension
//! 

use pyo3::exceptions;
use pyo3::prelude::PyResult;

use intervals::bounded::Bounded;

use intervals::nestedranges::NestedRanges;
/// Degrade a spatial coverage to a specific depth.
/// 
/// # Arguments
/// 
/// * ``coverage`` - The spatial coverage to degrade.
/// * ``depth`` - The depth to degrade the spatial coverage to.
/// 
/// # Errors
/// 
/// * ``depth`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 29]`
pub fn degrade_nested_ranges(coverage: &mut NestedRanges<u64>, depth: i8) -> PyResult<()> {
    if depth < 0 || depth > <u64>::MAXDEPTH {
        Err(exceptions::ValueError::py_err(
            format!("Depth must be in [0, {0}]", <u64>::MAXDEPTH)
        ))
    } else {
        coverage.degrade(depth);
        Ok(())
    }
}

use ndarray::{Axis, Array, Array1, Array2};
/// Flatten HEALPix cells to a specific depth
/// 
/// # Arguments
/// 
/// * ``coverage`` - The spatial coverage
/// * ``depth`` - The depth to flatten the coverage to.
pub fn flatten_pixels(data: Array2<u64>, depth: i8) -> Array1<u64> {
    let factor = 1 << (2 * (<u64>::MAXDEPTH - depth)) as u64;

    let flattened_intervals = &data / &Array::from_elem(data.shape(), factor);

    let mut flattened_pixels = Vec::<u64>::new();
    for interval in flattened_intervals.axis_iter(Axis(0)) {
        for pix in interval[0]..interval[1] {
            flattened_pixels.push(pix);
        }
    }
    Array1::from_vec(flattened_pixels)
}

use pyo3::types::{PyDict, PyString, PyList};
use pyo3::{ToPyObject, PyObject};
use pyo3::prelude::Python;

use std::ops::Range;
/// Deserialize a spatial coverage from a json python dictionary
/// 
/// JSON python dictionary stores (key, value) pair where:
/// 
/// * the ``key`` is a ``char`` being a depth
/// * the ``value`` is a list of HEALPix cell indices at the depth
///   indicated by the ``key``
/// 
/// # Arguments
/// 
/// * ``input`` - The input python dictionary
/// 
/// # Errors
/// 
/// * ``input`` dict must have string typed ``key``.
/// * ``input`` dict values must be a list of unsigned integer encoded
///   on 64 bits (i.e. an array of `u64`).
pub fn from_json(py: Python, input: &PyDict) -> PyResult<NestedRanges<u64>> {
    const TYPE_KEY_MSG_ERR: &'static str = "The key must be a python str that \n \
    encodes an integer on 1 byte. Ex: {'5': [0, 6, 7, ..., 9]}";
    const TYPE_VALUES_MSG_ERR: &'static str = "The values must be a list of unsigned \n \
    integer on 64 bits (8 bytes). Ex: {'5': [0, 6, 7, ..., 9]}";
    const EXTRACT_IPIX_FROM_LIST_MSG_ERR: &'static str = "Cannot extract 64 bits unsigned integer from the python list!";

    let mut ranges = Vec::<Range<u64>>::new();

    for (depth, pixels) in input.into_iter() {
        let depth = depth.downcast_ref::<PyString>()
            .map_err(|_| {
                exceptions::TypeError::py_err(TYPE_KEY_MSG_ERR)
            })?
            .to_string()?
            .parse::<i8>()
            .unwrap();

        let pixels = pixels.downcast_ref::<PyList>()
            .map_err(|_| {
                exceptions::TypeError::py_err(TYPE_VALUES_MSG_ERR)
            })?;

        let shift = ((<u64>::MAXDEPTH - depth) << 1) as u64;

        for p in pixels.into_iter() {
            let pixel = p
                .to_object(py)
                .extract::<u64>(py)
                .map_err(|_| {
                    exceptions::ValueError::py_err(EXTRACT_IPIX_FROM_LIST_MSG_ERR)
                })?;

            let e1 = pixel << shift;
            let e2 = (pixel + 1) << shift;
            ranges.push(e1..e2);
        }
    }

    let result = NestedRanges::<u64>::new(ranges)
        .make_consistent();
    Ok(result)
}

use std::collections::HashMap;
/// Serializes a spatial coverage to a JSON format
/// 
/// # Arguments
/// 
/// * ``coverage`` - The spatial coverage ranges to serialize.
pub fn to_json(py: Python, coverage: NestedRanges<u64>) -> PyResult<PyObject> {
    let result = PyDict::new(py);

    let size = (<u64>::MAXDEPTH + 1) as usize;
    let mut dict = HashMap::with_capacity(size);

    // Initialize the hash map with all the depths and
    // an empty vector of pixels associated to each of them
    for d in 0..size {
        dict.insert(d.to_string(), Vec::<u64>::new()).unwrap();
    }

    // Fill the hash map with the pixels
    for (depth, pix) in coverage.iter_depth_pix() {
        dict.get_mut(&depth.to_string())
            .unwrap()
            .push(pix);
    }

    // Fill the dictionary with only the depths
    // where pixels are found
    for (d, ipix) in &dict {
        if !ipix.is_empty() {
            result.set_item(d, ipix).map_err(|_| {
                exceptions::ValueError::py_err("An error occured when inserting items into the PyDict")
            })?;
        }
    }

    Ok(result.to_object(py))
}

/// Make a spatial coverage consistent
/// 
/// # Infos
/// 
/// This is an internal method whose purpose is not to be called
/// by an user. It is called inside of the `mocpy.IntervalSet` class.
/// 
/// # Arguments
/// 
/// * ``coverage`` - The spatial coverage to make consistent.
/// * ``min_depth`` - A minimum depth. This argument is optional.
///   A min depth means that we do not want any HEALPix cells to be
///   of depth < to ``min_depth``. This argument is used for example for
///   plotting a spatial coverage. All HEALPix cells of depth < 2 are splitted
///   into cells of depth 2.
/// 
/// # Errors
/// 
/// * ``min_depth`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 29]`
pub fn merge(coverage: NestedRanges<u64>, min_depth: i8) -> PyResult<NestedRanges<u64>> {
    // Convert the Array2<u64> to a NestedRanges<u64>
    // and make it consistent
    let mut coverage = coverage.make_consistent();

    // If a min depth has been given
    if min_depth != -1 {
        // Then we check its validity
        let max_depth = <u64>::MAXDEPTH;
        if min_depth < 0 || min_depth > max_depth {
            return Err(exceptions::ValueError::py_err("Min depth is not valid."));
        }

        // And perform the division of the ranges
        coverage = coverage.divide(min_depth);
    }

    Ok(coverage)
}

/// Compute the depth of a spatial coverage
/// 
/// # Arguments
/// 
/// * ``coverage`` - The input coverage.
pub fn depth(coverage: &NestedRanges<u64>) -> i8 {
    coverage.depth()
}

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