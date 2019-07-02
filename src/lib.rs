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

use ndarray::{Zip, Array, Array1, Array2, Axis};
use ndarray_parallel::prelude::*;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::{pymodule, Py, PyModule, PyResult, Python};
use pyo3::exceptions;
use pyo3::types::{PyDict, PyString, PyList};
use pyo3::{ToPyObject, PyObject};

use intervals::{NestedRanges, UniqRanges, RangesPy};
use intervals::bounded::Bounded;

use std::ops::Range;
use std::collections::HashMap;

mod ranges;

#[allow(unused_parens)]
#[pymodule]
fn core(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "from_lonlat")]
    fn from_lonlat_py(py: Python,
        depth: u8,
        lon: &PyArray1<f64>,
        lat: &PyArray1<f64>)
    -> Py<PyArray2<u64>> {
        let lon = lon.as_array()
            .to_owned()
            .into_raw_vec();
        let lat = lat.as_array()
            .to_owned()
            .into_raw_vec();

        let ranges = ranges::create_from_position(lon, lat, depth);
        
        let result: Array2<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "union")]
    fn union_py(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let a = a.as_array().to_owned();
        let b = b.as_array().to_owned();

        // Deal with the bound cases when
        // one of the two intervals are empty
        if a.is_empty() {
            return b
                .into_pyarray(py)
                .to_owned();
        } else if b.is_empty() {
            return a
                .into_pyarray(py)
                .to_owned();
        }

        let a: NestedRanges<u64> = RangesPy {
            data: a,
            min_depth: None,
            make_consistent: false,
        }.into();
        let b: NestedRanges<u64> = RangesPy {
            data: b,
            min_depth: None,
            make_consistent: false,
        }.into();

        let result = a.union(&b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "difference")]
    fn difference_py(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let a = a.as_array().to_owned();
        let b = b.as_array().to_owned();

        if a.is_empty() || b.is_empty() {
            return a.into_pyarray(py).to_owned();
        }

        let a: NestedRanges<u64> = RangesPy { 
            data: a,
            min_depth: None,
            make_consistent: false,
        }.into();
        let b: NestedRanges<u64> = RangesPy { 
            data: b,
            min_depth: None,
            make_consistent: false,
        }.into();

        let result = a.difference(&b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "intersection")]
    fn intersection_py(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let a = a.as_array().to_owned();
        let b = b.as_array().to_owned();

        if a.is_empty() {
            return a.into_pyarray(py).to_owned();
        } else if b.is_empty() {
            return b.into_pyarray(py).to_owned();
        }

        let a: NestedRanges<u64> = RangesPy { 
            data: a,
            min_depth: None,
            make_consistent: false,
        }.into();
        let b: NestedRanges<u64> = RangesPy { 
            data: b,
            min_depth: None,
            make_consistent: false,
        }.into();

        let result = a.intersection(&b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "complement")]
    fn complement_py(py: Python, a: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let a = a.as_array().to_owned();

        let ranges = if a.is_empty() {
            NestedRanges::<u64>::new(vec![], None, false)
        } else {
            RangesPy { 
                data: a,
                min_depth: None,
                make_consistent: false,
            }.into()
        };
        let result = ranges.complement();

        let result = if !result.is_empty() {
            result.into()
        } else {
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
        const EXTRACT_IPIX_FROM_LIST_MSG_ERR: &'static str = "Cannot extract 64bits unsigned integer from the python list!";

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

        let ranges = NestedRanges::<u64>::new(ranges, None, true);
        let result: Array2<u64> = ranges.into();
        Ok(result.into_pyarray(py).to_owned())
    }
    
    #[pyfn(m, "to_json")]
    fn to_json_py(py: Python, input: &PyArray2<u64>) -> PyResult<PyObject> {
        let input = input.as_array().to_owned();
        let result = PyDict::new(py);
        
        if !input.is_empty() {
            let ranges: NestedRanges<u64> = RangesPy { 
                data: input,
                min_depth: None,
                make_consistent: false,
            }.into();

            let size = (<u64>::MAXDEPTH + 1) as usize;
            let mut dict = HashMap::with_capacity(size);

            for d in 0..size {
                dict.insert(d.to_string(), Vec::<u64>::new()).unwrap();
            }

            for (depth, pix) in ranges.iter_depth_pix() {
                dict.get_mut(&depth.to_string())
                    .unwrap()
                    .push(pix);
            }

            for (d, ipix) in &dict {
                if !ipix.is_empty() {
                    result.set_item(d, ipix).map_err(|_| {
                        exceptions::ValueError::py_err("An error occured when inserting items into the PyDict")
                    })?;
                }
            }
        }

        Ok(result.to_object(py))
    }
    
    #[pyfn(m, "degrade")]
    fn degrade_py(py: Python, a: &PyArray2<u64>, depth: i8) -> Py<PyArray2<u64>> {
        let a = a.as_array().to_owned();
        if a.is_empty() {
            return a.into_pyarray(py).to_owned();
        }

        let mut ranges: NestedRanges<u64> = RangesPy {
            data: a,
            min_depth: None,
            make_consistent: true,
        }.into();

        ranges.degrade(depth);

        let result: Array2<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "merge_nested_intervals")]
    fn merge_nested_intervals_py(py: Python, a: &PyArray2<u64>, min_depth: i8) -> Py<PyArray2<u64>> {
        let a = a.as_array().to_owned();
        if a.is_empty() {
            return a.into_pyarray(py).to_owned();
        }
        let min_depth = if min_depth == -1 {
            None
        } else {
            Some(min_depth)
        };

        let ranges: NestedRanges<u64> = RangesPy {
            data: a,
            min_depth: min_depth,
            make_consistent: true,
        }.into();

        let result: Array2<u64> = ranges.into();
        result.to_owned().into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "depth")]
    fn depth_py(_py: Python, a: &PyArray2<u64>) -> PyResult<i8> {
        let a = a.as_array().to_owned();
        if a.is_empty() {
            return Ok(0);
        }
        let ranges: NestedRanges<u64> = RangesPy {
            data: a,
            min_depth: None,
            make_consistent: false,
        }.into();

        let depth = ranges.depth();
        Ok(depth)
    }
    
    #[pyfn(m, "to_nested")]
    fn to_nested_py(py: Python, uniq: &PyArray1<u64>) -> Py<PyArray2<u64>> {
        let uniq = uniq.as_array().to_owned();
        if uniq.is_empty() {
            return Array::zeros((1, 0)).into_pyarray(py).to_owned();
        }

        let shape = (uniq.shape()[0], 1);

        let start = uniq.into_shape(shape).unwrap();
        let end = &start + &Array::ones(shape);

        let input = stack![Axis(1), start, end].to_owned();

        let ranges: UniqRanges<u64> = RangesPy {
            data: input,
            min_depth: None,
            make_consistent: true,
        }.into();

        let ranges = ranges.to_nested();
        let result: Array2<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "to_uniq")]
    fn to_uniq_py(py: Python, nested: &PyArray2<u64>) -> Py<PyArray1<u64>> {
        let input = nested.as_array().to_owned();

        if input.is_empty() {
            return Array::zeros((0,)).into_pyarray(py).to_owned();
        }

        let ranges: NestedRanges<u64> = RangesPy {
            data: input,
            min_depth: None,
            make_consistent: false,
        }.into();

        let uniq_ranges = ranges.to_uniq();

        let mut data = Vec::<u64>::new();
        for r in uniq_ranges.iter() {
            for u in r.start..r.end {
                data.push(u);
            }
        }
        // Get a Array1 from the Vec<u64> without copying any data
        let result = Array1::from_vec(data);
        result.to_owned();

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
            let ranges: NestedRanges<u64> = RangesPy {
                data: ranges,
                min_depth: None,
                make_consistent: true,
            }.into();

            let result: Array2<u64> = ranges.into();
            Ok(result.into_pyarray(py).to_owned())
        }
    }
    
    #[pyfn(m, "flatten_pixels")]
    fn flatten_pixels_py(py: Python, nested: &PyArray2<u64>, depth: i8) -> Py<PyArray1<u64>> {
        let input = nested.as_array().to_owned();
        let factor = 1 << (2 * (<u64>::MAXDEPTH - depth)) as u64;

        let flattened_intervals = &input / &Array::from_elem(input.shape(), factor);

        let mut flattened_pixels = Vec::<u64>::new();
        for interval in flattened_intervals.axis_iter(Axis(0)) {
            for pix in interval[0]..interval[1] {
                flattened_pixels.push(pix);
            }
        }
        let result = Array1::from_vec(flattened_pixels).to_owned();
        result.into_pyarray(py).to_owned()
    }
    
    #[pyfn(m, "from_healpix_cells")]
    fn from_healpix_cells_py(py: Python, pixels: &PyArray1<u64>, depth: &PyArray1<i8>) -> Py<PyArray2<u64>> {
        let mut pixels = pixels.as_array().to_owned();
        let mut pixels_1 = &pixels + &Array::ones(pixels.shape());
        
        let depth = depth.as_array().to_owned();

        /*if pixels.shape() != depth.shape() {
            return Err(exceptions::IndexError::py_err("pixels and depth arrays must have the same shape"));
        }*/

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

        let pixel_ranges = stack![Axis(1), pixels, pixels_1].to_owned();

        let ranges: NestedRanges<u64> = RangesPy {
            data: pixel_ranges,
            min_depth: None,
            make_consistent: true,
        }.into();

        let result: Array2<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }
    
    Ok(())
}
