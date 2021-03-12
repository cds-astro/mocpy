#[cfg(feature = "rayon")]
extern crate intervals;
#[macro_use(concatenate)]

extern crate ndarray;
extern crate healpix;
extern crate num;
extern crate numpy;
extern crate rayon;
extern crate time;

extern crate pyo3;

#[macro_use]
extern crate lazy_static;

use ndarray::{Array, Array1, Array2, Axis};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::{pymodule, Py, PyModule, PyResult, Python};
use pyo3::types::{PyDict, PyList};
use pyo3::{PyObject, ToPyObject};

use intervals::nestedranges2d::NestedRanges2D;

use std::collections::HashMap;
use std::sync::Mutex;

pub mod coverage;
pub mod spatial_coverage;
pub mod temporal_coverage;
pub mod time_space_coverage;

type Coverage2DHashMap = HashMap<usize, NestedRanges2D<u64, u64>>;

lazy_static! {
    static ref COVERAGES_2D: Mutex<Coverage2DHashMap> = Mutex::new(HashMap::new());
    static ref NUM_COVERAGES_2D: Mutex<usize> = Mutex::new(0);
}

/// Insert a Time-Space coverage in the Hash Map
/// storing all the current 2D coverages
///
/// # Arguments
///
/// * `coverage` - The new Time-Space coverage to insert
///
/// # Panics
///
/// * This will panic if the `COVERAGES_2D` or `NUM_COVERAGES_2D`
///   are already held by the current thread
fn insert_new_coverage(coverage: NestedRanges2D<u64, u64>) -> usize {
    let mut coverages = COVERAGES_2D.lock().unwrap();
    let mut num_coverages = NUM_COVERAGES_2D.lock().unwrap();

    let index = *num_coverages;
    if let Some(_v) = coverages.insert(index, coverage) {
        panic!("There is already a coverage at this index.");
    }
    *num_coverages += 1;

    index
}

/// Remove a Time-Space coverage from the Hash Map
/// storing all the current 2D coverages
///
/// # Arguments
///
/// * `index` - The coverage to remove
///
/// # Panics
///
/// * If `COVERAGES_2D` is already held by the current thread.
fn remove_coverage(index: usize) {
    let mut coverages = COVERAGES_2D.lock().unwrap();
    let _coverage = coverages
        .remove(&index)
        // `None` is mapped to `Err(&'static str)`
        // because we suppose there should be a coverage
        // stored in the hash map at the `index` key.
        .expect("There is no coverage to remove");
}

/// Replace a Time-Space coverage at a specific index.
///
/// # Arguments
///
/// * `index` - The index of the Time-Space coverage to replace
/// * `coverage` - The new coverage
///
/// # Panics
///
/// * If no Time-Space coverage has been found in the hash map
/// for this specific `index`.
/// * If `COVERAGES_2D` is already held by the current thread.
fn update_coverage(index: usize, coverage: NestedRanges2D<u64, u64>) {
    let mut coverages = COVERAGES_2D.lock().unwrap();
    coverages
        .insert(index, coverage)
        // `None` is mapped to `Err(&'static str)`
        // because we suppose there should be a coverage
        // stored in the hash map at the `index` key.
        .expect("There is no coverage present");
}

#[pymodule]
fn mocpy(_py: Python, m: &PyModule) -> PyResult<()> {
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
    fn from_lonlat(
        py: Python,
        depth: i8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let lon = lon.as_array().to_owned().into_raw_vec();
        let lat = lat.as_array().to_owned().into_raw_vec();

        let ranges = spatial_coverage::create_from_position(lon, lat, depth)?;

        let result: Array2<u64> = ranges.into();
        Ok(result.into_pyarray(py).to_owned())
    }
    
    /// Create a 1D spatial coverage from a list of uniq cells each associated with a value.
    ///
    /// The coverage computed contains the cells summing from ``cumul_from`` to ``cumul_to``.
    ///
    /// # Arguments
    ///
    /// * ``uniq`` - Uniq HEALPix indices
    /// * ``values`` - Array containing the values associated for each cells.
    /// Must be of the same size of ``uniq`` and must sum to one.
    /// * ``cumul_from`` - The cumulative value from which cells are put in the coverage
    /// * ``cumul_to`` - The cumulative value to which cells are put in the coverage
    /// * ``max_depth`` -  the largest depth of the output MOC, which must be larger or equals to the largest
    ///   depth in the `uniq` values
    ///
    /// # Precondition
    ///
    /// * ``uniq`` and ``values`` must be of the same size
    /// * ``values`` must sum to one
    #[pyfn(m, "from_valued_hpx_cells")]
    fn from_valued_hpx_cells(
        py: Python,
        max_depth: u8,
        uniq: PyReadonlyArray1<u64>,
        values: PyReadonlyArray1<f64>,
        cumul_from: f64,
        cumul_to: f64,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let uniq = uniq.as_array().to_owned();
        let values = values.as_array().to_owned();

        let ranges = spatial_coverage::from_valued_healpix_cells(max_depth as u32, uniq, values, cumul_from, cumul_to)?;

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
    #[pyfn(m, "from_time_lonlat")]
    fn from_time_lonlat(
        index: usize,
        times: PyReadonlyArray1<f64>,
        d1: i8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
        d2: i8,
    ) -> PyResult<()> {
        let times = times.as_array()
            .to_owned()
            .into_raw_vec();
        let lon = lon.as_array()
            .to_owned()
            .into_raw_vec();
        let lat = lat.as_array()
            .to_owned()
            .into_raw_vec();

        let coverage = time_space_coverage::create_from_times_positions(times, lon, lat, d1, d2)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage);

        Ok(())
    }

    /// Create a 2D Time-Space coverage from a list of
    /// (time_range, longitude, latitude) tuples.
    ///
    /// # Arguments
    ///
    /// * ``times_min`` - The begining time of observation.
    /// * ``times_max`` - The ending time of observation.
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
    #[pyfn(m, "from_time_ranges_lonlat")]
    fn from_time_ranges_lonlat(
        index: usize,
        times_min: PyReadonlyArray1<f64>,
        times_max: PyReadonlyArray1<f64>,
        d1: i8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
        d2: i8,
    ) -> PyResult<()> {
        let times_min = times_min.as_array()
            .to_owned()
            .into_raw_vec();
        let times_max = times_max.as_array()
            .to_owned()
            .into_raw_vec();
        let lon = lon.as_array()
            .to_owned()
            .into_raw_vec();
        let lat = lat.as_array()
            .to_owned()
            .into_raw_vec();

        let coverage = time_space_coverage::create_from_time_ranges_positions(times_min, times_max, d1, lon, lat, d2)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage);

        Ok(())
    }

    /// Create a 2D Time-Space coverage from a list of
    /// (time_range, longitude, latitude, radius) tuples.
    ///
    /// # Arguments
    ///
    /// * ``times_min`` - The begining time of observation.
    /// * ``times_max`` - The ending time of observation.
    /// * ``d1`` - The depth along the Time axis.
    /// * ``lon`` - The longitudes in radians.
    /// * ``lat`` - The latitudes in radians.
    /// * ``radius`` - Radius in radians.
    /// * ``d2`` - The depth along the Space axis.
    ///
    /// # Precondition
    ///
    /// * ``lon``, ``lat`` and ``radius`` must be expressed in radians.
    /// * ``times`` must be expressed in jd.
    ///
    /// # Errors
    ///
    /// * ``lon``, ``lat``, ``times_min``, ``times_max`` and ``radius`` do not have the same length.
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 29]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    ///
    #[pyfn(m, "from_time_ranges_spatial_coverages")]
    fn from_time_ranges_spatial_coverages(
        py: Python,
        index: usize,
        times_min: PyReadonlyArray1<f64>,
        times_max: PyReadonlyArray1<f64>,
        d1: i8,
        spatial_coverages: &PyList,
    ) -> PyResult<()> {
        let times_min = times_min.as_array()
            .to_owned()
            .into_raw_vec();
        let times_max = times_max.as_array()
            .to_owned()
            .into_raw_vec();

        let coverage = time_space_coverage::from_time_ranges_spatial_coverages(py, times_min, times_max, d1, spatial_coverages)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage);

        Ok(())
    }

    #[pyfn(m, "project_on_first_dim")]
    fn project_on_first_dim(py: Python, ranges: PyReadonlyArray2<u64>, index: usize) -> Py<PyArray2<u64>> {
        // Build the input ranges from a Array2
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_nested_ranges_from_py(ranges);

        // Get the coverage and perform the projection
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            time_space_coverage::project_on_first_dim(&ranges, coverage)
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
    fn project_on_second_dim(
        py: Python,
        ranges: PyReadonlyArray2<u64>,
        index: usize,
    ) -> Py<PyArray2<u64>> {
        // Build the input ranges from a Array2
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_nested_ranges_from_py(ranges);

        // Get the coverage and perform the projection
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            time_space_coverage::project_on_second_dim(&ranges, coverage)
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
    fn coverage_2d_to_fits(py: Python, index: usize) -> Py<PyArray1<i64>> {
        // Get the coverage and flatten it
        // to a Array1
        let result: Array1<i64> = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            time_space_coverage::to_fits(coverage)
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
    /// This method returns a `PyValueError` if the `Array1` is not
    /// defined as above.
    #[pyfn(m, "coverage_2d_from_fits")]
    fn coverage_2d_from_fits(index: usize, data: PyReadonlyArray1<i64>) -> PyResult<()> {
        let data = data.as_array().to_owned();
        let coverage_from_fits = time_space_coverage::from_fits(data)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_fits);

        Ok(())
    }

    /// Create a new empty Time-Space coverage
    ///
    /// This method is called in the constructor of the
    /// `mocpy.STMOC` class
    ///
    /// # Returns
    ///
    /// The index of the newly created Time-Space coverage
    #[pyfn(m, "create_2d_coverage")]
    fn create_2d_coverage(_py: Python) -> usize {
        // Create new empty coverage
        let empty_coverage = time_space_coverage::new();

        // Insert a new coverage in the COVERAGES_2D
        // hash map and return its index key to python
        let result = insert_new_coverage(empty_coverage);

        result
    }

    /// Drop the content of a Time-Space coverage
    ///
    /// This method is automatically called by the
    /// Python garbage collector.
    #[pyfn(m, "drop_2d_coverage")]
    fn drop_2d_coverage(_py: Python, index: usize) {
        remove_coverage(index);
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
    fn coverage_2d_depth(_py: Python, index: usize) -> (i8, i8) {
        // Get the coverage and computes its depth
        // If the coverage is empty, the depth will be
        // (0, 0)
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            time_space_coverage::depth(coverage)
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
    fn coverage_2d_min_time(_py: Python, index: usize) -> PyResult<f64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        time_space_coverage::t_min(coverage)
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
    fn coverage_2d_max_time(_py: Python, index: usize) -> PyResult<f64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        time_space_coverage::t_max(coverage)
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
    fn coverage_2d_union(_py: Python, index: usize, id_left: usize, id_right: usize) {
        let result = {
            let coverages = COVERAGES_2D.lock().unwrap();
            // Get the left and right coverages
            let coverage_left = coverages.get(&id_left).unwrap();
            let coverage_right = coverages.get(&id_right).unwrap();

            // Perform the union
            let result = time_space_coverage::union(coverage_left, coverage_right);

            result
        };

        // Update the coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, result);
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
    fn coverage_2d_intersection(_py: Python, index: usize, id_left: usize, id_right: usize) {
        let result = {
            let coverages = COVERAGES_2D.lock().unwrap();
            // Get the left and right coverages
            let coverage_left = coverages.get(&id_left).unwrap();
            let coverage_right = coverages.get(&id_right).unwrap();

            // Perform the intersection
            let result = time_space_coverage::intersection(coverage_left, coverage_right);

            result
        };

        // Update the coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, result);
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
    fn coverage_2d_difference(_py: Python, index: usize, id_left: usize, id_right: usize) {
        let result = {
            let coverages = COVERAGES_2D.lock().unwrap();
            // Get the left and right coverages
            let coverage_left = coverages.get(&id_left).unwrap();
            let coverage_right = coverages.get(&id_right).unwrap();

            // Perform the difference
            let result = time_space_coverage::difference(coverage_left, coverage_right);

            result
        };

        // Update the coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, result)
    }

    /// Check the equality between two Time-Space coverages
    ///
    /// # Arguments
    ///
    /// * ``id_left`` - The index of the Time-Space coverage being
    ///   in the left of the operation.
    /// * ``id_right`` - The index of the Time-Space coverage being
    ///   in the right of the operation.
    #[pyfn(m, "coverage_2d_equality_check")]
    fn coverage_2d_equality_check(_py: Python, id_left: usize, id_right: usize) -> bool {
        let result = {
            let coverages = COVERAGES_2D.lock().unwrap();
            // Get the left and right coverages
            let cov_left = coverages.get(&id_left).unwrap();
            let cov_right = coverages.get(&id_right).unwrap();

            // Check the equality
            cov_left == cov_right
        };

        // Return the index of the newly created
        // coverage
        result
    }

    /// Checks whether a Time-Space coverage is empty.
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage to check
    ///   the emptiness.
    #[pyfn(m, "coverage_2d_is_empty")]
    fn coverage_2d_is_empty(_py: Python, index: usize) -> bool {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        coverage.is_empty()
    }

    /// Check if (time, position) tuples are contained into a Time-Space coverage
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage.
    /// * ``times`` - Times at which the positions have been observed.
    /// * ``lon`` - The longitudes.
    /// * ``lat`` - The latitudes.
    ///
    /// # Errors
    ///
    /// * If `lon`, `lat` and `times` do not have the same length
    #[pyfn(m, "coverage_2d_contains")]
    fn coverage_2d_contains(
        py: Python,
        index: usize,
        times: PyReadonlyArray1<f64>,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>) -> PyResult<Py<PyArray1<bool>>> {
        let times = times.as_array().to_owned();
        let lon = lon.as_array().to_owned();
        let lat = lat.as_array().to_owned();

        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        let mut result: Array1<bool> = Array::from_elem((lon.shape()[0],), false);
        time_space_coverage::contains(coverage, times, lon, lat, &mut result)?;
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Perform the union between two spatial coverages
    ///
    /// # Arguments
    ///
    /// * ``a`` - The spatial coverage being the left operand
    /// * ``b`` - The spatial coverage being the right operand
    #[pyfn(m, "coverage_union")]
    fn coverage_union(py: Python, a: PyReadonlyArray2<u64>, b: PyReadonlyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let cov_a = coverage::create_nested_ranges_from_py(ranges_a);
        let cov_b = coverage::create_nested_ranges_from_py(ranges_b);

        let result = spatial_coverage::union(&cov_a, &cov_b);

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
    fn coverage_difference(py: Python, a: PyReadonlyArray2<u64>, b: PyReadonlyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let cov_a = coverage::create_nested_ranges_from_py(ranges_a);
        let cov_b = coverage::create_nested_ranges_from_py(ranges_b);

        let result = spatial_coverage::difference(&cov_a, &cov_b);

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
    fn coverage_intersection(
        py: Python,
        a: PyReadonlyArray2<u64>,
        b: PyReadonlyArray2<u64>,
    ) -> Py<PyArray2<u64>> {
        let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let cov_a = coverage::create_nested_ranges_from_py(ranges_a);
        let cov_b = coverage::create_nested_ranges_from_py(ranges_b);

        let result = spatial_coverage::intersection(&cov_a, &cov_b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()
    }

    /// Computes the complement of the coverage
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The input spatial coverage
    #[pyfn(m, "coverage_complement")]
    fn coverage_complement(py: Python, ranges: PyReadonlyArray2<u64>) -> Py<PyArray2<u64>> {
        let ranges = ranges.as_array().to_owned();

        let coverage = coverage::create_nested_ranges_from_py(ranges);
        let result = spatial_coverage::complement(&coverage);

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
        let coverage = coverage::from_json(py, input)?;

        let result: Array2<u64> = coverage.into();
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Serialize a spatial coverage to a JSON format
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The spatial coverage ranges to serialize.
    #[pyfn(m, "coverage_to_json")]
    fn coverage_to_json(py: Python, ranges: PyReadonlyArray2<u64>) -> PyResult<PyObject> {
        let ranges = ranges.as_array().to_owned();

        let coverage = coverage::create_nested_ranges_from_py(ranges);

        let result = coverage::to_json(py, coverage)?;
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
    fn coverage_degrade(
        py: Python,
        ranges: PyReadonlyArray2<u64>,
        depth: i8,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = ranges.as_array().to_owned();

        let mut ranges = coverage::create_nested_ranges_from_py(ranges);
        coverage::degrade_nested_ranges(&mut ranges, depth)?;

        // Merge the overlapping intervals after degradation
        let ranges = ranges.make_consistent();

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
    fn coverage_merge_nested_intervals(
        py: Python,
        ranges: PyReadonlyArray2<u64>,
        min_depth: i8,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = ranges.as_array().to_owned();

        // Convert the Array2<u64> to a NestedRanges<u64>
        // and make it consistent
        let mut coverage = coverage::create_nested_ranges_from_py(ranges);
        coverage = coverage::merge(coverage, min_depth)?;

        let result: Array2<u64> = coverage.into();
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Compute the depth of a spatial coverage
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The input coverage.
    #[pyfn(m, "coverage_depth")]
    fn coverage_depth(_py: Python, ranges: PyReadonlyArray2<u64>) -> i8 {
        let ranges = ranges.as_array().to_owned();

        let coverage = coverage::create_nested_ranges_from_py(ranges);

        coverage::depth(&coverage)
    }

    /// Compute the sky fraction of a spatial coverage
    ///
    /// # Arguments
    ///
    /// * ``coverage`` - The spatial coverage
    /// * ``max_depth`` - The max depth of the spatial coverage.
    #[pyfn(m, "coverage_sky_fraction")]
    fn coverage_sky_fraction(_py: Python, ranges: PyReadonlyArray2<u64>) -> f32 {
        let ranges = ranges.as_array().to_owned();

        coverage::sky_fraction(&ranges)
    }

    /// Convert HEALPix cell indices from the **uniq** to the **nested** format.
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The HEALPix cells defined in the **uniq** format.
    #[pyfn(m, "to_nested")]
    fn to_nested(py: Python, ranges: PyReadonlyArray1<u64>) -> Py<PyArray2<u64>> {
        let ranges = ranges.as_array().to_owned();

        let result: Array2<u64> = if ranges.is_empty() {
            Array::zeros((1, 0))
        } else {
            let shape = (ranges.shape()[0], 1);

            let start = ranges.into_shape(shape).unwrap();
            let end = &start + &Array2::<u64>::ones(shape);

            let ranges = concatenate![Axis(1), start, end];
            let uniq_coverage = coverage::create_uniq_ranges_from_py(ranges).make_consistent();

            let nested_coverage = spatial_coverage::to_nested(uniq_coverage);
            nested_coverage.into()
        };

        result.into_pyarray(py).to_owned()
    }

    /// Convert HEALPix cell indices from the **nested** to the **uniq** format.
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The HEALPix cells defined in the **nested** format.
    #[pyfn(m, "to_uniq")]
    fn to_uniq(py: Python, ranges: PyReadonlyArray2<u64>) -> Py<PyArray1<u64>> {
        let ranges = ranges.as_array().to_owned();

        let result: Array1<u64> = if ranges.is_empty() {
            Array::zeros((0,))
        } else {
            let nested_coverage = coverage::create_nested_ranges_from_py(ranges);

            let uniq_coverage = nested_coverage.to_uniq();
            uniq_coverage.into()
        };

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
    fn from_time_ranges(
        py: Python,
        min_times: PyReadonlyArray1<f64>,
        max_times: PyReadonlyArray1<f64>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let min_times = min_times.as_array().to_owned();
        let max_times = max_times.as_array().to_owned();

        let coverage = temporal_coverage::from_time_ranges(min_times, max_times)?;

        let result: Array2<u64> = coverage.into();
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Flatten HEALPix cells to a specific depth
    ///
    /// # Arguments
    ///
    /// * ``data`` - The spatial coverage
    /// * ``depth`` - The depth to flatten the coverage to.
    #[pyfn(m, "flatten_pixels")]
    fn flatten_pixels(py: Python, data: PyReadonlyArray2<u64>, depth: i8) -> Py<PyArray1<u64>> {
        let data = data.as_array().to_owned();

        let result = coverage::flatten_pixels(data, depth);

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
    fn from_healpix_cells(
        py: Python,
        pixels: PyReadonlyArray1<u64>,
        depth: PyReadonlyArray1<i8>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let pixels = pixels.as_array().to_owned();
        let depth = depth.as_array().to_owned();

        let result = spatial_coverage::from_healpix_cells(pixels, depth)?;

        Ok(result.into_pyarray(py).to_owned())
    }

    Ok(())
}
