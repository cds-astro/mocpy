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

#[pymodule]
fn core(_py: Python, m: &PyModule) -> PyResult<()> {
    /// Create a 1D spatial coverage from a list of
    /// longitudes and latitudes
    /// 
    /// # Arguments
    /// 
    /// * ``depth`` - The depth of the coverage between `[0, <u64>::MAXDEPTH] = [0, 29]`
    /// * ``lon`` - The longitudes in radians
    /// * ``lat`` - The latitudes in radians
    /// 
    /// # Precondition
    /// 
    /// ``lon`` and ``lat`` must be expressed in radians and be valid.
    /// 
    /// # Errors
    /// 
    /// * ``lon`` and ``lat`` do not have the same length
    /// * ``depth`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 29]`
    #[pyfn(m, "from_lonlat")]
    fn from_lonlat_py(py: Python,
        depth: i8,
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

    /// Create a 2D Time-Space coverage from a list of
    /// (time, longitude, latitude) tuples.
    /// 
    /// # Arguments
    /// 
    /// * ``times`` - The times at which the sky coordinates have be given.
    /// * ``d1`` - The depth along the Time axis.
    /// * ``lon`` - The longitudes in radians
    /// * ``lat`` - The latitudes in radians
    /// * ``d2`` - The depth along the Space axis.
    /// 
    /// # Precondition
    /// 
    /// * ``lon`` and ``lat`` must be expressed in radians.
    /// * ``times`` must be expressed in jd.
    /// 
    /// # Errors
    /// 
    /// * ``lon``, ``lat`` and ``times`` do not have the same length.
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 29]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    /// 
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

    /// Project the Time-Space coverage into its first dimension
    /// (i.e. the Time axis)
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The constrained space set of ranges.
    /// * ``index`` - The index of the Time-Space coverage.
    /// 
    /// # Algorithm
    /// 
    /// Returns the union of the time ranges for which
    /// their space ranges is contained into ``y``.
    /// 
    /// # Panic
    /// 
    /// If the ``ranges`` is not valid i.e.:
    /// 
    /// * Contains ranges for which their inf bound is
    ///   superior to their sup bound.
    /// 
    /// This **should** not panic as this code is wrapped around MOCPy
    #[pyfn(m, "project_on_first_dim")]
    fn project_on_first_dim(py: Python,
        ranges: &PyArray2<u64>,
        index: usize,
    ) -> Py<PyArray2<u64>> {
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
        result.into_pyarray(py).to_owned()
    }

    /// Project the Time-Space coverage into its second dimension
    /// (i.e. the Space axis)
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The constrained time set of ranges.
    /// * ``index`` - The index of the Time-Space coverage.
    /// 
    /// # Algorithm
    /// 
    /// Returns the union of the spatial coverages for which
    /// their time ranges is contained into ``x``.
    /// 
    /// # Panic
    /// 
    /// If the ``ranges`` is not valid i.e.:
    /// 
    /// * Contains ranges for which their inf bound is
    ///   superior to their sup bound.
    /// 
    /// This **should** not panic as this code is wrapped around MOCPy
    #[pyfn(m, "project_on_second_dim")]
    fn project_on_second_dim(py: Python,
        ranges: &PyArray2<u64>,
        index: usize,
    ) -> Py<PyArray2<u64>> {
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
        result.into_pyarray(py).to_owned()
    }

    /// Serialize a Time-Space coverage into FITS
    /// 
    /// # Context
    /// 
    /// This is wrapped around the `serialize` method 
    /// of MOCPy to serialize a Time-Space coverage into
    /// FITS.
    /// 
    /// # Arguments 
    /// 
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize. 
    #[pyfn(m, "coverage_2d_to_fits")]
    fn coverage_2d_to_fits(py: Python,
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
    /// This method returns a `ValueError` if the `Array1` is not
    /// defined as above.
    use std::convert::TryFrom;
    #[pyfn(m, "coverage_2d_from_fits")]
    fn coverage_2d_from_fits(data: &PyArray1<i64>) -> PyResult<usize> {
        let data = data.as_array().to_owned();
        let new_coverage =  NestedRanges2D::<u64, u64>::try_from(data)
            .map_err(|msg| {
                    exceptions::ValueError::py_err(msg)
                })?;

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

    /// Computes the depth of a Time-Space coverage
    /// 
    /// # Arguments
    /// 
    /// * ``index`` - The index of the Time-Space coverage.
    /// 
    /// # Infos
    /// 
    /// If the Time-Space coverage is empty, the returned
    /// depth is `(0, 0)`.
    #[pyfn(m, "coverage_2d_depth")]
    fn coverage_2d_depth(_py: Python,
        index: usize,
    ) -> (i8, i8) {
        // Get the coverage and computes its depth
        // If the coverage is empty, the depth will be
        // (0, 0)
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            coverage.depth()
        };

        result
    }

    /// Returns the minimum time value of the Time-Space
    /// coverage
    /// 
    /// # Arguments
    /// 
    /// * ``index`` - The index of the Time-Space coverage.
    /// 
    /// # Errors
    /// 
    /// * If the coverage is empty.
    #[pyfn(m, "coverage_2d_min_time")]
    fn coverage_2d_min_time(_py: Python,
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

    /// Returns the maximum time value of the Time-Space
    /// coverage
    /// 
    /// # Arguments
    /// 
    /// * ``index`` - The index of the Time-Space coverage.
    /// 
    /// # Errors
    /// 
    /// * If the coverage is empty.
    #[pyfn(m, "coverage_2d_max_time")]
    fn coverage_2d_max_time(_py: Python,
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

    /// Perform the union between two Time-Space coverages.
    /// 
    /// # Arguments
    /// 
    /// * ``id_left`` - The index of the Time-Space coverage being
    ///   in the left of the operation.
    /// * ``id_right`` - The index of the Time-Space coverage being
    ///   in the right of the operation.
    #[pyfn(m, "coverage_2d_union")]
    fn coverage_2d_union(_py: Python,
        id_left: usize,
        id_right: usize,
    ) -> usize {
        let result = {
            let mut coverages = COVERAGES_2D.lock().unwrap();
            // Get the left and right coverages
            let cov_left = coverages.get(&id_left).unwrap();
            let cov_right = coverages.get(&id_right).unwrap();

            // Perform the union
            let result = cov_left.union(cov_right);

            // Insert the resulted coverage in the
            // hash map
            let index = coverages.len();
            coverages.insert(index, result);

            index
        };
        
        // Return the index of the newly created
        // coverage
        result
    }

    /// Perform the intersection between two Time-Space coverages.
    /// 
    /// # Arguments
    /// 
    /// * ``id_left`` - The index of the Time-Space coverage being
    ///   in the left of the operation.
    /// * ``id_right`` - The index of the Time-Space coverage being
    ///   in the right of the operation.
    #[pyfn(m, "coverage_2d_intersection")]
    fn coverage_2d_intersection(_py: Python,
        id_left: usize,
        id_right: usize,
    ) -> usize {
        let result = {
            let mut coverages = COVERAGES_2D.lock().unwrap();
            // Get the left and right coverages
            let cov_left = coverages.get(&id_left).unwrap();
            let cov_right = coverages.get(&id_right).unwrap();

            // Perform the intersection
            let result = cov_left.intersection(cov_right);

            // Insert the resulted coverage in the
            // hash map
            let index = coverages.len();
            coverages.insert(index, result);

            index
        };
        
        // Return the index of the newly created
        // coverage
        result
    }

    /// Perform the difference between two Time-Space coverages.
    /// 
    /// # Arguments
    /// 
    /// * ``id_left`` - The index of the Time-Space coverage being
    ///   in the left of the operation.
    /// * ``id_right`` - The index of the Time-Space coverage being
    ///   in the right of the operation.
    #[pyfn(m, "coverage_2d_difference")]
    fn coverage_2d_difference(_py: Python,
        id_left: usize,
        id_right: usize,
    ) -> usize {
        let result = {
            let mut coverages = COVERAGES_2D.lock().unwrap();
            // Get the left and right coverages
            let cov_left = coverages.get(&id_left).unwrap();
            let cov_right = coverages.get(&id_right).unwrap();

            // Perform the intersection
            let result = cov_left.difference(cov_right);

            // Insert the resulted coverage in the
            // hash map
            let index = coverages.len();
            coverages.insert(index, result);

            index
        };
        
        // Return the index of the newly created
        // coverage
        result
    }

    /// Perform the union between two spatial coverages
    /// 
    /// # Arguments
    /// 
    /// * ``a`` - The spatial coverage being the left operand
    /// * ``b`` - The spatial coverage being the right operand
    #[pyfn(m, "coverage_union")]
    fn coverage_union(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let ranges_a = ranges::create_nested_ranges_from_py(ranges_a);
        let ranges_b = ranges::create_nested_ranges_from_py(ranges_b);

        let result = ranges_a.union(&ranges_b);

        let result: Array2<u64> = result.into();
        result.to_owned().into_pyarray(py).to_owned()
    }

    /// Perform the difference between two spatial coverages
    /// 
    /// # Arguments
    /// 
    /// * ``a`` - The spatial coverage being the left operand
    /// * ``b`` - The spatial coverage being the right operand
    #[pyfn(m, "coverage_difference")]
    fn coverage_difference(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let ranges_a = ranges::create_nested_ranges_from_py(ranges_a);
        let ranges_b = ranges::create_nested_ranges_from_py(ranges_b);

        let result = ranges_a.difference(&ranges_b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }

    /// Perform the intersection between two spatial coverages
    /// 
    /// # Arguments
    /// 
    /// * ``a`` - The spatial coverage being the left operand
    /// * ``b`` - The spatial coverage being the right operand
    #[pyfn(m, "coverage_intersection")]
    fn coverage_intersection(py: Python, a: &PyArray2<u64>, b: &PyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let ranges_a = ranges::create_nested_ranges_from_py(ranges_a);
        let ranges_b = ranges::create_nested_ranges_from_py(ranges_b);

        let result = ranges_a.intersection(&ranges_b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }
    
    /// Computes the complement of the coverage
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The input spatial coverage
    #[pyfn(m, "coverage_complement")]
    fn coverage_complement(py: Python, ranges: &PyArray2<u64>) -> Py<PyArray2<u64>> {
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
    #[pyfn(m, "coverage_from_json")]
    fn coverage_from_json(py: Python, input: &PyDict) -> PyResult<Py<PyArray2<u64>>> {
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
    
    /// Serialize a spatial coverage to a JSON format
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The spatial coverage ranges to serialize.
    #[pyfn(m, "coverage_to_json")]
    fn coverage_to_json(py: Python, ranges: &PyArray2<u64>) -> PyResult<PyObject> {
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
    #[pyfn(m, "coverage_degrade")]
    fn coverage_degrade(py: Python, ranges: &PyArray2<u64>, depth: i8) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = ranges.as_array().to_owned();

        let mut ranges = ranges::create_nested_ranges_from_py(ranges);
        ranges::degrade_nested_ranges(&mut ranges, depth)?;

        let result: Array2<u64> = ranges.into();
        Ok(result.into_pyarray(py).to_owned())
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
    /// * ``ranges`` - The spatial coverage ranges to make consistent.
    /// * ``min_depth`` - A minimum depth. This argument is optional.
    ///   A min depth means that we do not want any HEALPix cells to be
    ///   of depth < to ``min_depth``. This argument is used for example for
    ///   plotting a spatial coverage. All HEALPix cells of depth < 2 are splitted
    ///   into cells of depth 2.
    /// 
    /// # Errors
    /// 
    /// * ``min_depth`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 29]`
    #[pyfn(m, "coverage_merge_nested_intervals")]
    fn coverage_merge_nested_intervals(py: Python, ranges: &PyArray2<u64>, min_depth: i8) -> PyResult<Py<PyArray2<u64>>> {
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

    /// Compute the depth of a spatial coverage
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The input coverage.
    #[pyfn(m, "coverage_depth")]
    fn coverage_depth(_py: Python, ranges: &PyArray2<u64>) -> i8 {
        let ranges = ranges.as_array().to_owned();

        let ranges = ranges::create_nested_ranges_from_py(ranges);

        let depth = ranges.depth();
        depth
    }

    /// Convert HEALPix cell indices from the **uniq** to the **nested** format.
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The HEALPix cells defined in the **uniq** format.
    #[pyfn(m, "to_nested")]
    fn to_nested(py: Python, ranges: &PyArray1<u64>) -> Py<PyArray2<u64>> {
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

    /// Convert HEALPix cell indices from the **nested** to the **uniq** format.
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The HEALPix cells defined in the **nested** format.
    #[pyfn(m, "to_uniq")]
    fn to_uniq(py: Python, ranges: &PyArray2<u64>) -> Py<PyArray1<u64>> {
        let ranges = ranges.as_array().to_owned();

        if ranges.is_empty() {
            return Array::zeros((0,)).into_pyarray(py).to_owned();
        }

        let ranges = ranges::create_nested_ranges_from_py(ranges);

        let ranges = ranges.to_uniq();
        let result: Array1<u64> = ranges.into();
        result.into_pyarray(py).to_owned()
    }

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

    /// Flatten HEALPix cells to a specific depth
    /// 
    /// # Arguments
    /// 
    /// * ``ranges`` - The spatial coverage
    /// * ``depth`` - The depth to flatten the coverage to.
    #[pyfn(m, "flatten_pixels")]
    fn flatten_pixels(py: Python, nested: &PyArray2<u64>, depth: i8) -> Py<PyArray1<u64>> {
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

    /// Create a spatial coverage from a list of HEALPix cell indices.
    /// 
    /// # Arguments
    /// 
    /// * ``pixels`` - A set of HEALPix cell indices
    /// * ``depth`` - The depths of each HEALPix cell indices
    /// 
    /// # Precondition
    /// 
    /// ``pixels`` and ``depth`` must be valid. This means that:
    /// 
    /// * ``depth`` contains values in the range `[0, <T>::MAXDEPTH] = [0, 29]`
    /// * ``pixels`` contains values in the range `[0, 12*4**(depth)]`
    /// 
    /// # Errors
    /// 
    /// * ``depth`` and ``pixels`` have not the same length.
    #[pyfn(m, "from_healpix_cells")]
    fn from_healpix_cells(py: Python, pixels: &PyArray1<u64>, depth: &PyArray1<i8>) -> PyResult<Py<PyArray2<u64>>> {
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
