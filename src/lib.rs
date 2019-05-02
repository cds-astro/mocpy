extern crate intervals;

extern crate ndarray;
extern crate ndarray_parallel;
extern crate numpy;
extern crate healpix;
extern crate num;
extern crate time;

extern crate pyo3;

use ndarray::{Zip, ArrayD, ArrayViewD, ArrayViewMutD, Array1, Array2, ArrayViewMut2};
use ndarray_parallel::prelude::*;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayDyn};
use pyo3::prelude::{pymodule, Py, PyModule, PyResult, Python};
use num::{Integer, PrimInt};

use intervals::intervals::Intervals;
use intervals::ranges::Ranges;
use intervals::bounded::Bounded;

use time::PreciseTime;

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
        
        let start = PreciseTime::now();
        // pix is full we can compute the merged intervals
        let intervals = {
            let ptr = pixels.as_mut_ptr() as *mut Range<u64>;
            let cap = len;

            mem::forget(pixels);

            let data = unsafe {
                Vec::from_raw_parts(ptr, len, cap)
            };

            data.into()
        };
        let end = PreciseTime::now();
        println!("{} seconds for allocation", start.to(end));

        let result = intervals_to_2darray(intervals);
        result.into_pyarray(py).to_owned()
    }

    /// Suppose the input array is contiguous in memory
    unsafe fn array2d_to_intervals(mut input: Array2<u64>) -> Intervals<u64> {
        let len = input.shape()[0];
        let cap = len;
        let ptr = input.as_mut_ptr() as *mut Range<u64>;

        mem::forget(input);

        let ranges = Vec::from_raw_parts(ptr, len, cap);

        Intervals::Nested(
            Ranges::<u64>::new(ranges)
        )
    }

    fn intervals_to_2darray<T: Integer + PrimInt + Bounded<T> + Send + 'static>(intervals: Intervals<T>) -> Array2<T> {
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
    /*
    #[pyfn(m, "union")]
    fn union_py(py: Python, input1: &PyArray2<u64>, input2: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let a = input1.as_array();
        let b = input2.as_array();

        let intervals1 = unsafe {
            ndarray_to_intervals::<u64>(a)
        };
        let intervals2 = unsafe {
            ndarray_to_intervals::<u64>(b)
        };

        intervals1.union(intervals2);

        a.into_pyarray().to_owned()
    }
    */
    #[pyfn(m, "depth")]
    fn depth_py(py: Python, input: &PyArray2<u64>) -> PyResult<i8> {
        let input = input.as_array().to_owned();
        let intervals = unsafe {
            array2d_to_intervals(input)
        };
        
        let depth = intervals.depth();
        Ok(depth)
    }

    Ok(())
}