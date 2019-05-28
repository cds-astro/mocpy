extern crate intervals;
#[macro_use(stack)]
extern crate ndarray;
extern crate ndarray_parallel;
extern crate numpy;
extern crate healpix;
extern crate num;
extern crate time;

extern crate pyo3;

use ndarray::{Zip, Array, Array1, Array2, Axis};
use ndarray_parallel::prelude::*;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::{pymodule, Py, PyModule, PyResult, Python};
use pyo3::exceptions;
use pyo3::types::{PyDict, PyString, PyList};
use pyo3::ToPyObject;
use num::{Integer, PrimInt};

use intervals::intervals::Intervals;
use intervals::ranges::Ranges;
use intervals::bounded::Bounded;

//use time::PreciseTime;

use std::ops::Range;
use std::mem;

#[pymodule]
fn core(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "from_lonlat")]
    fn from_lonlat_py(py: Python,
        depth: u8,
        lon: &PyArray1<f64>,
        lat: &PyArray1<f64>)
    -> Py<PyArray2<u64>> {
        let lon = lon.as_array();
        let lat = lat.as_array();

        let len = lon.shape()[0];

        // Allocate a new result array that will
        // store the merged nested intervals of the MOC
        // No copy of the data occuring
        let ranges: Vec<Range<u64>> = vec![0..1; len];
        let mut pixels = Array1::from_vec(ranges);

        let shift = (u64::MAXDEPTH as u8 - depth) << 1;

        let layer = healpix::nested::get_or_create(depth);
        Zip::from(&mut pixels)
            .and(&lon)
            .and(&lat)
            .par_apply(|p, &lon, &lat| {
                let pix = layer.hash(lon, lat);

                let e1 = pix << shift;
                let e2 = (pix + 1) << shift;
                *p = e1..e2;
            });

        //let start = PreciseTime::now();
        let intervals = {
            let ptr = pixels.as_mut_ptr() as *mut Range<u64>;
            let cap = len;

            mem::forget(pixels);

            let data = unsafe {
                Vec::from_raw_parts(ptr, len, cap)
            };

            data.into()
        };
        //let end = PreciseTime::now();
        //println!("{} time elapsed", start.to(end));

        let result = intervals_to_2darray(intervals);
        result.into_pyarray(py).to_owned()
    }

    /// Suppose the input array is contiguous in memory
    unsafe fn array2d_to_intervals(mut input: Array2<u64>, nested: bool, min_depth: Option<i8>) -> Intervals<u64> {
        let len = input.shape()[0];
        let cap = len;
        let ptr = input.as_mut_ptr() as *mut Range<u64>;

        mem::forget(input);

        let ranges = Vec::from_raw_parts(ptr, len, cap);

        if nested {
            Intervals::Nested(
                Ranges::<u64>::new(ranges, min_depth)
            )
        } else {
            Intervals::Uniq(
                Ranges::<u64>::new(ranges, min_depth)
            )
        }
    }

    fn intervals_to_2darray<T>(intervals: Intervals<T>) -> Array2<T>
    where T : Integer + PrimInt + Bounded<T> + Send + std::fmt::Debug + 'static {
        if let Intervals::Nested(ranges) = intervals {
            // Cast Vec<Range<u64>> to Vec<u64>
            let len = ranges.0.len();
            let data = ranges.to_flat_vec();

            // Get a Array1 from the Vec<u64> without copying any data
            let result = Array1::from_vec(data);

            // Reshape the result to get a Array2 of shape (N x 2) where N is the number 
            // of HEALPix cell contained in the moc
            result.into_shape((len, 2))
                  .unwrap()
                  .to_owned()
        } else {
            unreachable!()
        }
    }
    
    #[pyfn(m, "union")]
    fn union_py(py: Python, input1: &PyArray2<u64>, input2: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let input1 = input1.as_array().to_owned();
        let input2 = input2.as_array().to_owned();

        let mut i1 = unsafe {
            array2d_to_intervals(input1, true, None)
        };
        let i2 = unsafe {
            array2d_to_intervals(input2, true, None)
        };

        i1.union(i2);

        let result = intervals_to_2darray(i1);
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "difference")]
    fn difference_py(py: Python, input1: &PyArray2<u64>, input2: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let input1 = input1.as_array().to_owned();
        let input2 = input2.as_array().to_owned();

        let mut i1 = unsafe {
            array2d_to_intervals(input1, true, None)
        };
        let i2 = unsafe {
            array2d_to_intervals(input2, true, None)
        };

        i1.difference(i2);

        let result = intervals_to_2darray(i1);
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "intersection")]
    fn intersection_py(py: Python, input1: &PyArray2<u64>, input2: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let input1 = input1.as_array().to_owned();
        let input2 = input2.as_array().to_owned();

        let mut i1 = unsafe {
            array2d_to_intervals(input1, true, None)
        };
        let i2 = unsafe {
            array2d_to_intervals(input2, true, None)
        };

        i1.intersection(i2);

        let result = intervals_to_2darray(i1);
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "complement")]
    fn complement_py(py: Python, input: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let input = input
            .as_array()
            .to_owned();

        let mut intervals = unsafe {
            array2d_to_intervals(input, true, None)
        };
        intervals.complement();

        let result = intervals_to_2darray(intervals);
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

        let intervals = Intervals::Nested(
            Ranges::<u64>::new(ranges, None)
        );
        let result = intervals_to_2darray(intervals);
        Ok(result.into_pyarray(py).to_owned())
    }

    #[pyfn(m, "degrade")]
    fn degrade_py(py: Python, input: &PyArray2<u64>, depth: i8) -> Py<PyArray2<u64>> {
        let input = input.as_array().to_owned();

        let mut intervals = unsafe {
            array2d_to_intervals(input, true, None)
        };

        intervals.degrade(depth);

        let result = intervals_to_2darray(intervals);
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "refine_to_order")]
    fn refine_to_order_py(py: Python, input: &PyArray2<u64>, min_depth: i8) -> Py<PyArray2<u64>> {
        let input = input.as_array().to_owned();

        let intervals = unsafe {
            array2d_to_intervals(input, true, Some(min_depth))
        };

        let result = intervals_to_2darray(intervals);
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "depth")]
    fn depth_py(_py: Python, input: &PyArray2<u64>) -> PyResult<i8> {
        let input = input.as_array().to_owned();
        let intervals = unsafe {
            array2d_to_intervals(input, true, None)
        };
        
        let depth = intervals.depth();
        Ok(depth)
    }

    #[pyfn(m, "to_nested")]
    fn to_nested_py(py: Python, uniq: &PyArray1<u64>) -> Py<PyArray2<u64>> {
        let start = uniq.as_array().to_owned();
        let shape = (start.shape()[0], 1);

        let start = start.into_shape(shape).unwrap();
        let end = &start + &Array::ones(shape);

        let input = stack![Axis(1), start, end].to_owned();

        let intervals = unsafe {
            array2d_to_intervals(input, false, None)
        };

        let intervals = intervals.to_nested();
        let result = intervals_to_2darray(intervals);
        result.into_pyarray(py).to_owned()
    }

    #[pyfn(m, "to_uniq")]
    fn to_uniq_py(py: Python, nested: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let input = nested.as_array().to_owned();

        let intervals = unsafe {
            array2d_to_intervals(input, true, None)
        };

        let intervals = intervals.to_uniq();

        let result = intervals_to_2darray(intervals);
        result.into_pyarray(py).to_owned()
    }

    Ok(())
}
