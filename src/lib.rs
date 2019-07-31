extern crate intervals;
#[macro_use(stack)]

extern crate ndarray;
extern crate ndarray_parallel;
extern crate numpy;
extern crate healpix;
extern crate num;
extern crate time;
extern crate rayon;

extern crate pyo3;

#[macro_use]
extern crate lazy_static;

use ndarray::{Zip, Array, Array1, Array2, Axis};
use ndarray_parallel::prelude::*;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::{pymodule, Py, PyModule, PyResult, Python};
use pyo3::exceptions;
use pyo3::types::{PyDict, PyString, PyList};
use pyo3::{ToPyObject, PyObject};

use intervals::nestedranges::NestedRanges;

use intervals::nestedranges2d::NestedRanges2D;
use intervals::bounded::Bounded;

use std::ops::Range;
use std::collections::HashMap;
use std::sync::Mutex;

mod ranges;

type Coverage2DHashMap = HashMap<usize, NestedRanges2D<u64, u64>>;

lazy_static! {
    static ref COVERAGES_2D: Mutex<Coverage2DHashMap> = Mutex::new(HashMap::new());
}

#[allow(unused_parens)]
#[pymodule]
fn core(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "from_lonlat")]
    fn from_lonlat_py(py: Python,
        depth: u8,
        lon: &PyArray1<f64>,
        lat: &PyArray1<f64>)
    -> PyResult<Py<PyArray2<u64>>> {
        let lon = lon.as_array()
            .to_owned()
            .into_raw_vec();
        let lat = lat.as_array()
            .to_owned()
            .into_raw_vec();

        let ranges = ranges::create_from_position(lon, lat, depth)?;
        let result: Array2<u64> = ranges.into();
        Ok(result.into_pyarray(py).to_owned())
    }

    #[pyfn(m, "from_time_lonlat")]
    fn from_time_lonlat_py(
        times: &PyArray1<f64>,
        d1: i8,
        lon: &PyArray1<f64>,
        lat: &PyArray1<f64>,
        d2: i8)
    -> PyResult<usize> {
        let times = times.as_array()
            .to_owned()
            .into_raw_vec();
        let lon = lon.as_array()
            .to_owned()
            .into_raw_vec();
        let lat = lat.as_array()
            .to_owned()
            .into_raw_vec();

        let ts_ranges = ranges::create_from_time_position(times, lon, lat, d1, d2)?;

        // Insert a new coverage in the COVERAGES_2D
        // hash map and return its index key to python
        let result = {
            let mut d = COVERAGES_2D.lock().unwrap();
            
            let index = d.len();
            d.insert(index, ts_ranges);
            
            index
        };

        Ok(result)
    }

    #[pyfn(m, "project_on_first_dim")]
    fn project_on_first_dim(py: Python,
        ranges: &PyArray2<u64>,
        index: usize,
    ) -> PyResult<Py<PyArray2<u64>>> {
        // Build the input ranges from a Array2
        let ranges = ranges.as_array().to_owned();
        let ranges = ranges::create_nested_ranges_from_py(ranges);

        // Get the coverage and perform the projection
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            NestedRanges2D::project_on_first_dim(&ranges, coverage)
        };

        // Convert the result back to an ndarray::Array2
        let result: Array2<u64> = result.into();
        Ok(result.into_pyarray(py).to_owned())
    }

    #[pyfn(m, "project_on_second_dim")]
    fn project_on_second_dim(py: Python,
        ranges: &PyArray2<u64>,
        index: usize,
    ) -> PyResult<Py<PyArray2<u64>>> {
        // Build the input ranges from a Array2
        let ranges = ranges.as_array().to_owned();
        let ranges = ranges::create_nested_ranges_from_py(ranges);
        // Get the coverage and perform the projection
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            NestedRanges2D::project_on_second_dim(&ranges, coverage)
        };

        // Convert the result back to an ndarray::Array2
        let result: Array2<u64> = result.into();
        Ok(result.into_pyarray(py).to_owned())
    }

    #[pyfn(m, "ranges2d_to_fits")]
    fn ranges2d_to_fits(py: Python,
        index: usize,
    ) -> Py<PyArray1<i64>> {
        // Get the coverage and flatten it
        // to a Array1
        let result: Array1<i64> = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            coverage.into()
        };

        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "ranges2d_from_fits")]
    fn ranges2d_from_fits(data: &PyArray1<i64>) -> PyResult<usize> {
        let data = data.as_array().to_owned();
        let new_coverage: NestedRanges2D<u64, u64> = data.into();
        
        // Insert a new coverage in the COVERAGES_2D
        // hash map and return its index key to python
        let result = {
            let mut d = COVERAGES_2D.lock().unwrap();
            
            let index = d.len();
            d.insert(index, new_coverage);
            
            index
        };

        Ok(result)
    }

    #[pyfn(m, "ranges2d_depth")]
    fn ranges2d_depth(_py: Python,
        index: usize,
    ) -> PyResult<(i8, i8)> {
        // Get the coverage and computes its depth
        // If the coverage is empty, the depth will be
        // (0, 0)
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            coverage.depth()
        };

        Ok(result)
    }

    #[pyfn(m, "ranges2d_min_time")]
    fn ranges2d_min_time(_py: Python,
        index: usize,
    ) -> PyResult<f64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        let t_min = coverage.t_min()
            .map_err(|msg| {
                exceptions::ValueError::py_err(msg)
            })?;

        Ok((t_min as f64) / 86400000000_f64)
    }

    #[pyfn(m, "ranges2d_max_time")]
    fn ranges2d_max_time(_py: Python,
        index: usize,
    ) -> PyResult<f64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        let t_max = coverage.t_max()
            .map_err(|msg| {
                exceptions::ValueError::py_err(msg)
            })?;
        
        Ok((t_max as f64) / 86400000000_f64)
    }

    #[pyfn(m, "union")]
    fn union_py(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let ranges_a = ranges::create_nested_ranges_from_py(ranges_a);
        let ranges_b = ranges::create_nested_ranges_from_py(ranges_b);

        let result = ranges_a.union(&ranges_b);

        let result: Array2<u64> = result.into();
        result.to_owned().into_pyarray(py).to_owned()
    }

    #[pyfn(m, "difference")]
    fn difference_py(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let ranges_a = ranges::create_nested_ranges_from_py(ranges_a);
        let ranges_b = ranges::create_nested_ranges_from_py(ranges_b);

        let result = ranges_a.difference(&ranges_b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "intersection")]
    fn intersection_py(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let ranges_a = ranges::create_nested_ranges_from_py(ranges_a);
        let ranges_b = ranges::create_nested_ranges_from_py(ranges_b);

        let result = ranges_a.intersection(&ranges_b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "complement")]
    fn complement_py(py: Python, ranges: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges = ranges.as_array().to_owned();

        let ranges = ranges::create_nested_ranges_from_py(ranges);
        let result = ranges.complement();

        let result = if !result.is_empty() {
            result.into()
        } else {
            // TODO: try without this condition
            Array::zeros((1, 0))
        };

        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "from_json")]
    fn from_json_py(py: Python, input: &PyDict) -> PyResult<Py<PyArray2<u64>>> {
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

        let ranges = NestedRanges::<u64>::new(ranges)
            .make_consistent();
        let result: Array2<u64> = ranges.into();
        Ok(result.into_pyarray(py).to_owned())
    }
    
    #[pyfn(m, "to_json")]
    fn to_json_py(py: Python, ranges: &PyArray2<u64>) -> PyResult<PyObject> {
        let ranges = ranges.as_array().to_owned();
        let result = PyDict::new(py);
        
        let ranges = ranges::create_nested_ranges_from_py(ranges);

        let size = (<u64>::MAXDEPTH + 1) as usize;
        let mut dict = HashMap::with_capacity(size);

        // Initialize the hash map with all the depths and
        // an empty vector of pixels associated to each of them
        for d in 0..size {
            dict.insert(d.to_string(), Vec::<u64>::new()).unwrap();
        }

        // Fill the hash map with the pixels
        for (depth, pix) in ranges.iter_depth_pix() {
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

    #[pyfn(m, "degrade")]
    fn degrade_py(py: Python, ranges: &PyArray2<u64>, depth: i8) -> Py<PyArray2<u64>> {
        let ranges = ranges.as_array().to_owned();

        let mut ranges = ranges::create_nested_ranges_from_py(ranges)
            .make_consistent();
        ranges.degrade(depth);

        let result: Array2<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "merge_nested_intervals")]
    fn merge_nested_intervals_py(py: Python, ranges: &PyArray2<u64>, min_depth: i8) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = ranges.as_array().to_owned();

        // Convert the Array2<u64> to a NestedRanges<u64>
        // and make it consistent
        let mut ranges = ranges::create_nested_ranges_from_py(ranges)
            .make_consistent();
        
        // If a min depth has been given
        if min_depth != -1 {
            // Then we check its validity
            let max_depth = <u64>::MAXDEPTH;
            if min_depth < 0 || min_depth > max_depth {
                return Err(exceptions::ValueError::py_err("Min depth is not valid."));
            }

            // And perform the division of the ranges
            ranges = ranges.divide(min_depth);
        }

        let result: Array2<u64> = ranges.into();
        Ok(result.into_pyarray(py).to_owned())
    }

    #[pyfn(m, "depth")]
    fn depth_py(_py: Python, ranges: &PyArray2<u64>) -> i8 {
        let ranges = ranges.as_array().to_owned();

        let ranges = ranges::create_nested_ranges_from_py(ranges);

        let depth = ranges.depth();
        depth
    }

    #[pyfn(m, "to_nested")]
    fn to_nested_py(py: Python, ranges: &PyArray1<u64>) -> Py<PyArray2<u64>> {
        let ranges = ranges.as_array().to_owned();
        if ranges.is_empty() {
            return Array::zeros((1, 0)).into_pyarray(py).to_owned();
        }

        let shape = (ranges.shape()[0], 1);

        let start = ranges.into_shape(shape).unwrap();
        let end = &start + &Array::ones(shape);

        let ranges = stack![Axis(1), start, end];
        let ranges = ranges::create_uniq_ranges_from_py(ranges)
            .make_consistent();

        let ranges = ranges.to_nested();
        let result: Array2<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "to_uniq")]
    fn to_uniq_py(py: Python, ranges: &PyArray2<u64>) -> Py<PyArray1<u64>> {
        let ranges = ranges.as_array().to_owned();

        if ranges.is_empty() {
            return Array::zeros((0,)).into_pyarray(py).to_owned();
        }

        let ranges = ranges::create_nested_ranges_from_py(ranges);

        let ranges = ranges.to_uniq();
        let result: Array1<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "from_time_ranges")]
    fn from_time_ranges(py: Python, min_times: &PyArray1<f64>, max_times: &PyArray1<f64>) -> PyResult<Py<PyArray2<u64>>> {
        let min_times = min_times.as_array();
        let max_times = max_times.as_array();

        if min_times.shape() != max_times.shape() {
            Err(exceptions::ValueError::py_err("min and max ranges have not the same shape"))
        } else {
            if min_times.is_empty() {
                return Ok(Array::zeros((1, 0))
                    .into_pyarray(py)
                    .to_owned());
            }
            let shape = (min_times.shape()[0], 1);

            let min_times = min_times
                .into_shape(shape)
                .unwrap();
            let max_times = max_times
                .into_shape(shape)
                .unwrap();

            let min_times = &min_times * &Array::from_elem(shape, 86400000000_f64);
            let min_times = min_times.mapv(|e| e as u64);

            let max_times = &max_times * &Array::from_elem(shape, 86400000000_f64);
            let max_times = max_times.mapv(|e| e as u64 + 1);

            let ranges = stack![Axis(1), min_times, max_times].to_owned();
            let ranges = ranges::create_nested_ranges_from_py(ranges)
                .make_consistent();

            let result: Array2<u64> = ranges.into();
            Ok(result.into_pyarray(py).to_owned())
        }
    }

    #[pyfn(m, "flatten_pixels")]
    fn flatten_pixels_py(py: Python, nested: &PyArray2<u64>, depth: i8) -> Py<PyArray1<u64>> {
        let input = nested.as_array();
        let factor = 1 << (2 * (<u64>::MAXDEPTH - depth)) as u64;

        let flattened_intervals = &input / &Array::from_elem(input.shape(), factor);

        let mut flattened_pixels = Vec::<u64>::new();
        for interval in flattened_intervals.axis_iter(Axis(0)) {
            for pix in interval[0]..interval[1] {
                flattened_pixels.push(pix);
            }
        }
        let result = Array1::from_vec(flattened_pixels);
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "from_healpix_cells")]
    fn from_healpix_cells_py(py: Python, pixels: &PyArray1<u64>, depth: &PyArray1<i8>) -> PyResult<Py<PyArray2<u64>>> {
        let mut pixels = pixels.as_array().to_owned();
        let mut pixels_1 = &pixels + &Array::ones(pixels.shape());
        
        let depth = depth.as_array();

        if pixels.shape() != depth.shape() {
            return Err(exceptions::IndexError::py_err("pixels and depth arrays must have the same shape"));
        }

        Zip::from(&mut pixels)
            .and(&mut pixels_1)
            .and(&depth)
            .par_apply(|pix, pix1, &d| {
                let factor = 2 * (<u64>::MAXDEPTH - d);
                *pix <<= factor;
                *pix1 <<= factor;
            });

        let shape = (pixels.shape()[0], 1);
        let pixels = pixels.into_shape(shape).unwrap();
        let pixels_1 = pixels_1.into_shape(shape).unwrap();

        let ranges = stack![Axis(1), pixels, pixels_1].to_owned();

        let ranges = ranges::create_nested_ranges_from_py(ranges)
            .make_consistent();

        let result: Array2<u64> = ranges.into();
        Ok(result.into_pyarray(py).to_owned())
    }
    
    Ok(())
}
