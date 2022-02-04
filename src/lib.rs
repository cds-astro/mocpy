#[cfg(feature = "rayon")]
extern crate intervals;

extern crate ndarray;
extern crate healpix;
extern crate num;
extern crate numpy;
extern crate rayon;
extern crate time;

extern crate pyo3;

#[macro_use]
extern crate lazy_static;

use std::ops::Range;
use std::sync::Mutex;
use std::path::Path;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use ndarray::{Array, Array1, Array2, Ix2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};

use pyo3::prelude::{pymodule, Py, PyModule, PyResult, Python};
use pyo3::types::{PyDict, PyList, PyString};
use pyo3::{PyObject, ToPyObject, exceptions};

use moc::qty::{MocQty, Hpx};
use moc::elemset::range::{
    MocRanges,
    uniq::HpxUniqRanges
};
use moc::moc::{CellMOCIterator, CellMOCIntoIterator};
use moc::moc::range::RangeMOC;
use moc::ranges::{SNORanges, Ranges};
use moc::hpxranges2d::TimeSpaceMoc;

pub mod ndarray_fromto;
pub mod coverage;
pub mod spatial_coverage;
pub mod temporal_coverage;
pub mod time_space_coverage;

use crate::ndarray_fromto::{ranges_to_array2, mocranges_to_array2, vec_range_to_array2}; // uniq_ranges_to_array1

type Coverage2DHashMap = HashMap<usize, TimeSpaceMoc<u64, u64>>;

lazy_static! {
    static ref COVERAGES_2D: Mutex<Coverage2DHashMap> = Mutex::new(HashMap::new());
    static ref NUM_COVERAGES_2D: AtomicUsize = AtomicUsize::new(0);
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
fn insert_new_coverage(coverage: TimeSpaceMoc<u64, u64>) -> usize {
    let mut coverages = COVERAGES_2D.lock().unwrap();
    let index = NUM_COVERAGES_2D.fetch_add(1, Ordering::SeqCst );
    if let Some(_v) = coverages.insert(index, coverage) {
        panic!("There is already a coverage at this index.");
    }
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
fn update_coverage(index: usize, coverage: TimeSpaceMoc<u64, u64>) {
    let mut coverages = COVERAGES_2D.lock().unwrap();
    coverages
        .insert(index, coverage)
        // `None` is mapped to `Err(&'static str)`
        // because we suppose there should be a coverage
        // stored in the hash map at the `index` key.
        .expect("There is no coverage present");
}

fn coverage_op<O>(py: Python, a: PyReadonlyArray2<u64>, b: PyReadonlyArray2<u64>, op: O)
                  -> Py<PyArray2<u64>>
    where
      O: Fn(Ranges<u64>, Ranges<u64>) -> Ranges<u64>
{
    let ranges_a = a.as_array().to_owned();
    let ranges_b = b.as_array().to_owned();

    let cov_a = coverage::create_ranges_from_py_unchecked(ranges_a);
    let cov_b = coverage::create_ranges_from_py_unchecked(ranges_b);

    let result = op(cov_a, cov_b);

    let result: Array2<u64> = ranges_to_array2(result);
    result.to_owned().into_pyarray(py).to_owned()
}

fn coverage_complement<Q, F>(py: Python, ranges: PyReadonlyArray2<u64>, to_moc_ranges: F) -> Py<PyArray2<u64>>
    where
      Q: MocQty<u64>,
      F: Fn(Array2<u64>) -> MocRanges<u64, Q>
{
    let ranges = ranges.as_array().to_owned();

    let coverage = to_moc_ranges(ranges);
    let result = coverage.complement();

    let result = mocranges_to_array2(result);

    result.into_pyarray(py).to_owned()
}


fn coverage_degrade<Q, F>(
    py: Python,
    ranges: PyReadonlyArray2<u64>,
    depth: u8,
    to_moc_ranges: F,
) -> PyResult<Py<PyArray2<u64>>>
    where
      Q: MocQty<u64>,
      F: Fn(Array2<u64>) -> MocRanges<u64, Q>
{
    let ranges = ranges.as_array().to_owned();
    let mut ranges = to_moc_ranges(ranges);
    coverage::degrade_ranges(&mut ranges, depth)?;
    // The result is already consistent
    let result: Array2<u64> = mocranges_to_array2(ranges);
    Ok(result.into_pyarray(py).to_owned())
}

fn coverage_merge_intervals<Q, F>(
    py: Python,
    ranges: PyReadonlyArray2<u64>,
    min_depth: i8,
    to_moc_ranges: F,
) -> PyResult<Py<PyArray2<u64>>>
    where
      Q: MocQty<u64>,
      F: Fn(Array2<u64>) -> MocRanges<u64, Q>
{
    let ranges = ranges.as_array().to_owned();
    let mut coverage = to_moc_ranges(ranges);
    coverage = coverage::merge(coverage, min_depth)?;
    let result: Array2<u64> = mocranges_to_array2(coverage);
    Ok(result.into_pyarray(py).to_owned())
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
        depth: u8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let lon = lon.as_array().to_owned().into_raw_vec();
        let lat = lat.as_array().to_owned().into_raw_vec();

        let ranges = spatial_coverage::create_from_position(lon, lat, depth)?;

        let result: Array2<u64> = mocranges_to_array2(ranges);
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
    /// * `asc`: cumulative value computed from lower to highest densities instead of from highest to lowest
    /// * `strict`: (sub-)cells overlapping the `cumul_from` or `cumul_to` values are not added
    /// * `no_split`: cells overlapping the `cumul_from` or `cumul_to` values are not recursively split
    /// * `reverse_decent`: perform the recursive decent from the highest cell number to the lowest (to be compatible with Aladin)
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
        asc: bool,
        strict: bool,
        no_split: bool,
        reverse_decent: bool,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let uniq = uniq.as_array().to_owned();
        let values = values.as_array().to_owned();
        
        let ranges = spatial_coverage::from_valued_healpix_cells_with_opt(
            max_depth,
            uniq,
            values,
            cumul_from,
            cumul_to,
            asc,
            strict,
            no_split,
            reverse_decent,
        )?; //from_valued_healpix_cells(max_depth, uniq, values, cumul_from, cumul_to)?;

        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Create a 2D Time-Space coverage from a list of
    /// (time, longitude, latitude) tuples.
    ///
    /// # Arguments
    ///
    /// * ``times`` - The times at which the sky coordinates have be given, in jd coded
    ///   on doubles (=> not precise to the microsecond).
    /// * ``d1`` - The depth along the Time axis.
    /// * ``lon`` - The longitudes in radians
    /// * ``lat`` - The latitudes in radians
    /// * ``d2`` - The depth along the Space axis.
    ///
    /// # Precondition
    ///
    /// * ``lon`` and ``lat`` must be expressed in radians.
    /// * ``times`` must be expressed in jd, on doubles (=> not precise to the microsecond).
    ///
    /// # Errors
    ///
    /// * ``lon``, ``lat`` and ``times`` do not have the same length.
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 61]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    ///
    /// # Remark
    ///
    /// Method kept temporarily to ensure backward compatibility.
    ///
    #[pyfn(m, "from_time_lonlat_approx")]
    fn from_time_lonlat_approx(
        index: usize,
        times: PyReadonlyArray1<f64>,
        d1: u8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
        d2: u8,
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

        let coverage = time_space_coverage::create_from_times_positions_approx(times, lon, lat, d1, d2)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage);

        Ok(())
    }

    /// Create a 2D Time-Space coverage from a list of
    /// (time, longitude, latitude) tuples.
    ///
    /// # Arguments
    ///
    /// * ``times`` - The times at which the sky coordinates have be given, in microsecond since JD=0.
    /// * ``d1`` - The depth along the Time axis.
    /// * ``lon`` - The longitudes in radians
    /// * ``lat`` - The latitudes in radians
    /// * ``d2`` - The depth along the Space axis.
    ///
    /// # Precondition
    ///
    /// * ``lon`` and ``lat`` must be expressed in radians.
    /// * ``times`` must be expressed in jd, in microsecond since JD=0.
    ///
    /// # Errors
    ///
    /// * ``lon``, ``lat`` and ``times`` do not have the same length.
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 61]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    #[pyfn(m, "from_time_lonlat")]
    fn from_time_lonlat(
        index: usize,
        times: PyReadonlyArray1<u64>,
        d1: u8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
        d2: u8,
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
    /// * ``times_min`` - The begining time of observation, in jd coded
    ///                   on doubles (=> not precise to the microsecond).
    /// * ``times_max`` - The ending time of observation, in jd coded
    ///                   on doubles (=> not precise to the microsecond).
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
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 61]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    ///
    /// # Remark
    ///
    /// Method kept temporarily to ensure backward compatibility.
    ///
    #[pyfn(m, "from_time_ranges_lonlat_approx")]
    fn from_time_ranges_lonlat_approx(
        index: usize,
        times_min: PyReadonlyArray1<f64>,
        times_max: PyReadonlyArray1<f64>,
        d1: u8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
        d2: u8,
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

        let coverage = time_space_coverage::create_from_time_ranges_positions_approx(times_min, times_max, d1, lon, lat, d2)?;

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
    /// * ``times_min`` - The begining time of observation, in microsecond since JD=0.
    /// * ``times_max`` - The ending time of observation, in microsecond since JD=0.
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
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 61]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    ///
    #[pyfn(m, "from_time_ranges_lonlat")]
    fn from_time_ranges_lonlat(
        index: usize,
        times_min: PyReadonlyArray1<u64>,
        times_max: PyReadonlyArray1<u64>,
        d1: u8,
        lon: PyReadonlyArray1<f64>,
        lat: PyReadonlyArray1<f64>,
        d2: u8,
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
    /// * ``times_min`` - The begining time of observation in jd coded
    ///                   on doubles (=> not precise to the microsecond).
    /// * ``times_max`` - The ending time of observation, in jd coded
    ///                   on doubles (=> not precise to the microsecond).
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
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 61]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    ///
    /// # Remark
    ///
    /// Method kept temporarily to ensure backward compatibility.
    ///
    #[pyfn(m, "from_time_ranges_spatial_coverages_approx")]
    fn from_time_ranges_spatial_coverages_approx(
        py: Python,
        index: usize,
        times_min: PyReadonlyArray1<f64>,
        times_max: PyReadonlyArray1<f64>,
        d1: u8,
        spatial_coverages: &PyList,
    ) -> PyResult<()> {
        let times_min = times_min.as_array()
            .to_owned()
            .into_raw_vec();
        let times_max = times_max.as_array()
            .to_owned()
            .into_raw_vec();

        let coverage = time_space_coverage::from_time_ranges_spatial_coverages_approx(py, times_min, times_max, d1, spatial_coverages)?;

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
    /// * ``times_min`` - The begining time of observation, in microsecond since JD=0.
    /// * ``times_max`` - The ending time of observation, in microsecond since JD=0.
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
    /// * ``d1`` is not comprised in `[0, <T>::MAXDEPTH] = [0, 61]`
    /// * ``d2`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
    ///
    #[pyfn(m, "from_time_ranges_spatial_coverages")]
    fn from_time_ranges_spatial_coverages(
        py: Python,
        index: usize,
        times_min: PyReadonlyArray1<u64>,
        times_max: PyReadonlyArray1<u64>,
        d1: u8,
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
        let ranges = coverage::create_hpx_ranges_from_py_unchecked(ranges);

        // Get the coverage and perform the projection
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            time_space_coverage::project_on_first_dim(&ranges, coverage)
        };

        // Convert the result back to an ndarray::Array2
        let result: Array2<u64> = mocranges_to_array2(result);
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
        let ranges = coverage::create_time_ranges_from_py_uncheked(ranges);

        // Get the coverage and perform the projection
        let result = {
            let res = COVERAGES_2D.lock().unwrap();
            let coverage = res.get(&index).unwrap();

            time_space_coverage::project_on_second_dim(&ranges, coverage)
        };

        // Convert the result back to an ndarray::Array2
        let result: Array2<u64> = mocranges_to_array2(result);
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

    /// Serialize a Time-Space coverage into a FITS file
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize.
    /// * ``path`` - the path of the output file
    #[pyfn(m, "coverage_2d_to_fits_file")]
    fn coverage_2d_to_fits_file(path: String, index: usize) -> PyResult<()> {
        // Get the coverage and flatten it
        // to a Array1
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();
        let(depth_max_t, depth_max_s) = coverage.compute_min_depth();
        time_space_coverage::to_fits_file(depth_max_t, depth_max_s, coverage, Path::new(&path))
          .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Serialize a Time-Space coverage into an ASCII file
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize.
    /// * ``path`` - the path of the output file
    #[pyfn(m, "coverage_2d_to_ascii_file")]
    fn coverage_2d_to_ascii_file(path: String, index: usize) -> PyResult<()> {
        // Get the coverage and flatten it
        // to a Array1
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();
        let(depth_max_t, depth_max_s) = coverage.compute_min_depth();
        time_space_coverage::to_ascii_file(depth_max_t, depth_max_s, coverage, path)
          .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Serialize a Time-Space coverage into an ASCII string
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize.
    /// * ``path`` - the path of the output file
    #[pyfn(m, "coverage_2d_to_ascii_str")]
    fn coverage_2d_to_ascii_str(py: Python, index: usize) -> Py<PyString> {
        // Get the coverage and flatten it
        // to a Array1
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();
        let(depth_max_t, depth_max_s) = coverage.compute_min_depth();
        PyString::new(py, &time_space_coverage::to_ascii_str(depth_max_t, depth_max_s, coverage)).into()
    }

    /// Serialize a Time-Space coverage into a JSON file
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize.
    /// * ``path`` - the path of the output file
    #[pyfn(m, "coverage_2d_to_json_file")]
    fn coverage_2d_to_json_file(path: String, index: usize) -> PyResult<()> {
        // Get the coverage and flatten it
        // to a Array1
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();
        let(depth_max_t, depth_max_s) = coverage.compute_min_depth();
        time_space_coverage::to_json_file(depth_max_t, depth_max_s, coverage, path)
          .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Serialize a Time-Space coverage into a JSON file
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize.
    /// * ``path`` - the path of the output file
    #[pyfn(m, "coverage_2d_to_json_str")]
    fn coverage_2d_to_json_str(py: Python, index: usize) -> Py<PyString> {
        // Get the coverage and flatten it
        // to a Array1
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();
        let(depth_max_t, depth_max_s) = coverage.compute_min_depth();
        PyString::new(py, &time_space_coverage::to_json_str(depth_max_t, depth_max_s, coverage)).into()
    }

    /*
    /// Deserialize a Time-Space coverage from a JSON file
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize.
    /// * ``path`` - the path of the output file
    #[pyfn(m, "coverage_2d_from_json_file")]
    fn coverage_2d_from_json_file(path: String, index: usize) -> PyResult<()> {
        // Get the coverage and flatten it
        // to a Array1
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();
        let(depth_max_t, depth_max_s) = coverage.compute_min_depth();
        time_space_coverage::to_json_file(depth_max_t, depth_max_s, coverage, path)
          .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Deserialize a Time-Space coverage into a JSON file
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage
    ///   to serialize.
    /// * ``path`` - the path of the output file
    #[pyfn(m, "coverage_2d_from_json_str")]
    fn coverage_2d_from_json_str(py: Python, index: usize) -> Py<PyString> {
        // Get the coverage and flatten it
        // to a Array1
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();
        let(depth_max_t, depth_max_s) = coverage.compute_min_depth();
        PyString::new(py, &time_space_coverage::to_json_str(depth_max_t, depth_max_s, coverage)).into()
    }*/




    /// Deserialize a Time-Space coverage from FITS using the pre v2.0 MOC standard.
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
    #[pyfn(m, "coverage_2d_from_fits_pre_v2")]
    fn coverage_2d_from_fits_pre_v2(index: usize, data: PyReadonlyArray1<i64>) -> PyResult<()> {
        let data = data.as_array().to_owned();
        let coverage_from_fits = time_space_coverage::from_fits_pre_v2(data)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_fits);

        Ok(())
    }

    /// Deserialize a Time-Space coverage from FITS using the v2.0 MOC standard.
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
    /// Its memory layout contains a list of time ranges followed by the
    /// list of space ranges referred to that time ranges.
    /// The most significant bit (MSB) of time ranges bounds is set to one so that one can
    /// distinguish them from space ranges.
    /// This is different from a negative value because we do not use the two's complement
    /// representation, only a flag set on the MSB.
    ///
    /// This method returns a `PyValueError` if the `Array1` is not
    /// defined as above.
    #[pyfn(m, "coverage_2d_from_fits")]
    fn coverage_2d_from_fits(index: usize, data: PyReadonlyArray1<u64>) -> PyResult<()> {
        let data = data.as_array().to_owned();
        let coverage_from_fits = time_space_coverage::from_fits(data)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_fits);

        Ok(())
    }

    /// Deserialize a Time-Space coverage from a FITS file (compatible with the MOC v2.0 standard).
    ///
    /// # Arguments
    ///
    /// * ``index`` - the index used to store the Time-Space coverage
    /// * ``path`` - the FITS file path
    ///
    /// # Warning
    ///
    /// This function is not compatible with pre-v2.0 MOC standard.
    ///
    /// # Errors
    ///
    /// This method returns a `PyIOError` if the the function fails in writing the FITS file.
    #[pyfn(m, "coverage_2d_from_fits_file")]
    fn coverage_2d_from_fits_file(index: usize, path: String) -> PyResult<()> {
        let coverage_from_fits = time_space_coverage::from_fits_file(Path::new(&path))?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_fits);

        Ok(())
    }

    /// Deserialize a Time-Space coverage from an ASCII file (compatible with the MOC v2.0 standard).
    ///
    /// # Arguments
    ///
    /// * ``index`` - the index used to store the Time-Space coverage
    /// * ``path`` - the ASCII file path
    ///
    /// # Errors
    ///
    /// This method returns a `PyIOError` if the the function fails in writing the FITS file.
    #[pyfn(m, "coverage_2d_from_ascii_file")]
    fn coverage_2d_from_ascii_file(index: usize, path: String) -> PyResult<()> {
        let coverage_from_ascii = time_space_coverage::from_ascii_file(Path::new(&path))?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_ascii);

        Ok(())
    }

    /// Deserialize a Time-Space coverage from an JSON file.
    ///
    /// # Arguments
    ///
    /// * ``index`` - the index used to store the Time-Space coverage
    /// * ``path`` - the JSON file path
    ///
    /// # Errors
    ///
    /// This method returns a `PyIOError` if the the function fails in writing the FITS file.
    #[pyfn(m, "coverage_2d_from_json_file")]
    fn coverage_2d_from_json_file(index: usize, path: String) -> PyResult<()> {
        let coverage_from_json = time_space_coverage::from_json_file(Path::new(&path))?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_json);

        Ok(())
    }

    /// Deserialize a Time-Space coverage from an ASCII string (compatible with the MOC v2.0 standard).
    ///
    /// # Arguments
    ///
    /// * ``index`` - the index used to store the Time-Space coverage
    /// * ``ascii`` - the ASCII string
    ///
    /// # Errors
    ///
    /// This method returns a `PyIOError` if the the function fails in writing the FITS file.
    #[pyfn(m, "coverage_2d_from_ascii_str")]
    fn coverage_2d_from_ascii_str(index: usize, ascii: String) -> PyResult<()> {
        let coverage_from_ascii = time_space_coverage::from_ascii_str(ascii)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_ascii);

        Ok(())
    }

    /// Deserialize a Time-Space coverage from an JSON string.
    ///
    /// # Arguments
    ///
    /// * ``index`` - the index used to store the Time-Space coverage
    /// * ``json`` - the JSON string
    ///
    /// # Errors
    ///
    /// This method returns a `PyIOError` if the the function fails in writing the FITS file.
    #[pyfn(m, "coverage_2d_from_json_str")]
    fn coverage_2d_from_json_str(index: usize, json: String) -> PyResult<()> {
        let coverage_from_json = time_space_coverage::from_json_str(json)?;

        // Update a coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, coverage_from_json);

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
        insert_new_coverage(empty_coverage)
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
    fn coverage_2d_depth(_py: Python, index: usize) -> (u8, u8) {
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
    ///
    /// # Remark
    ///
    /// Method kept temporarily to ensure backward compatibility.
    ///
    #[pyfn(m, "coverage_2d_min_time_approx")]
    fn coverage_2d_min_time_approx(_py: Python, index: usize) -> PyResult<f64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        time_space_coverage::t_min_jd(coverage)
    }

    /// Returns the minimum time value of the Time-Space
    /// coverage, in microarcsec since jd=0.
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage.
    ///
    /// # Errors
    ///
    /// * If the coverage is empty.
    ///
    #[pyfn(m, "coverage_2d_min_time")]
    fn coverage_2d_min_time(_py: Python, index: usize) -> PyResult<u64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        time_space_coverage::t_min_mircosecond_since_jd_org(coverage)
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
    ///
    /// # Remark
    ///
    /// Method kept temporarily to ensure backward compatibility.
    ///
    #[pyfn(m, "coverage_2d_max_time_approx")]
    fn coverage_2d_max_time_approx(_py: Python, index: usize) -> PyResult<f64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        time_space_coverage::t_max_jd(coverage)
    }

    /// Returns the maximum time value of the Time-Space
    /// coverage, in microseconds since jd=0.
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage.
    ///
    /// # Errors
    ///
    /// * If the coverage is empty.
    ///
    #[pyfn(m, "coverage_2d_max_time")]
    fn coverage_2d_max_time(_py: Python, index: usize) -> PyResult<u64> {
        // Get the coverage
        let res = COVERAGES_2D.lock().unwrap();
        let coverage = res.get(&index).unwrap();

        time_space_coverage::t_max_mircosecond_since_jd_org(coverage)
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
            time_space_coverage::union(coverage_left, coverage_right)
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
            time_space_coverage::intersection(coverage_left, coverage_right)
        };

        // Update the coverage in the COVERAGES_2D
        // hash map and return its index key to python
        update_coverage(index, result)
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
            time_space_coverage::difference(coverage_left, coverage_right)
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
    ///
    /// # Remark
    ///
    /// Method kept temporarily to ensure backward compatibility.
    ///
    #[pyfn(m, "coverage_2d_contains_approx")]
    fn coverage_2d_contains_approx(
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
        time_space_coverage::contains_approx(coverage, times, lon, lat, &mut result)?;
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Check if (time, position) tuples are contained into a Time-Space coverage
    ///
    /// # Arguments
    ///
    /// * ``index`` - The index of the Time-Space coverage.
    /// * ``times`` - Times at which the positions have been observed, in microsec since jd=0
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
        times: PyReadonlyArray1<u64>,
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

    /// Perform the union between two generic coverages
    ///
    /// # Arguments
    ///
    /// * ``a`` - The spatial coverage being the left operand
    /// * ``b`` - The spatial coverage being the right operand
    #[pyfn(m, "coverage_union")]
    fn coverage_union(py: Python, a: PyReadonlyArray2<u64>, b: PyReadonlyArray2<u64>) -> Py<PyArray2<u64>> {
        /*let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let cov_a = coverage::create_ranges_from_py_unchecked(ranges_a);
        let cov_b = coverage::create_ranges_from_py_unchecked(ranges_b);

        let result = cov_a.union(&cov_b);

        let result: Array2<u64> = result.into();
        result.to_owned().into_pyarray(py).to_owned()*/
        coverage_op(py, a, b, |cov_a, cov_b| cov_a.union(&cov_b))
    }



    /// Perform the difference between two generic coverages
    ///
    /// # Arguments
    ///
    /// * ``a`` - The spatial coverage being the left operand
    /// * ``b`` - The spatial coverage being the right operand
    #[pyfn(m, "coverage_difference")]
    fn coverage_difference(py: Python, a: PyReadonlyArray2<u64>, b: PyReadonlyArray2<u64>) -> Py<PyArray2<u64>> {
        /*let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let cov_a = coverage::create_ranges_from_py(ranges_a);
        let cov_b = coverage::create_ranges_from_py(ranges_b);

        let result = cov_a.difference(&cov_b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()*/
        coverage_op(py, a, b, |cov_a, cov_b| cov_a.difference(&cov_b))
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
        /*let ranges_a = a.as_array().to_owned();
        let ranges_b = b.as_array().to_owned();

        let cov_a = coverage::create_ranges_from_py(ranges_a);
        let cov_b = coverage::create_ranges_from_py(ranges_b);

        let result = cov_a.intersection(&cov_b);

        let result: Array2<u64> = result.into();
        result.into_pyarray(py).to_owned()*/
        coverage_op(py, a, b, |cov_a, cov_b| cov_a.intersection(&cov_b))
    }

    /// Computes the complement of the given nested/ring coverage
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The input spatial coverage
    #[pyfn(m, "hpx_coverage_complement")]
    fn hpx_coverage_complement(py: Python, ranges: PyReadonlyArray2<u64>) -> Py<PyArray2<u64>> {
        coverage_complement(py, ranges, coverage::create_hpx_ranges_from_py_unchecked)
    }

    /// Computes the complement of the given time coverage
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The input time coverage
    #[pyfn(m, "time_coverage_complement")]
    fn time_coverage_complement(py: Python, ranges: PyReadonlyArray2<u64>) -> Py<PyArray2<u64>> {
        coverage_complement(py, ranges, coverage::create_time_ranges_from_py_uncheked)
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

        let result: Array2<u64> = mocranges_to_array2(coverage);
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

        let coverage = coverage::create_hpx_ranges_from_py_unchecked(ranges);

        let result = coverage::to_json(py, coverage)?;
        Ok(result.to_object(py))
    }


    /// Serialize a spatial MOC into an FITS file.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    /// * ``path`` - The file path
    #[pyfn(m, "spatial_moc_to_fits_file")]
    fn spatial_moc_to_fits_file(
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
        path: String,
    ) -> PyResult<()> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        spatial_coverage::to_fits_file(depth, ranges, path)
          .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Serialize a spatial MOC into an ASCII file.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    /// * ``path`` - The file path
    #[pyfn(m, "spatial_moc_to_ascii_file")]
    fn spatial_moc_to_ascii_file(
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
        path: String,
    ) -> PyResult<()> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        spatial_coverage::to_ascii_file(depth, ranges, path)
          .map_err(exceptions::PyIOError::new_err)
    }

    /// Serialize a spatial MOC into a ASCII string.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    #[pyfn(m, "spatial_moc_to_ascii_str")]
    fn spatial_moc_to_ascii_str(
        py: Python,
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
    ) -> Py<PyString> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        PyString::new(py,&spatial_coverage::to_ascii_str(depth, ranges)).into()
    }

    /// Serialize a spatial MOC into a JSON file.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    /// * ``path`` - The file path
    #[pyfn(m, "spatial_moc_to_json_file")]
    fn spatial_moc_to_json_file(
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
        path: String,
    ) -> PyResult<()> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        spatial_coverage::to_json_file(depth, ranges, path)
          .map_err(exceptions::PyIOError::new_err)
    }

    /// Serialize a spatial MOC into a JSON string.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    #[pyfn(m, "spatial_moc_to_json_str")]
    fn spatial_moc_to_json_str(
        py: Python,
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
    ) -> Py<PyString> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        PyString::new(py,&spatial_coverage::to_json_str(depth, ranges)).into()
    }

    /// Create a 1D spatial coverage from the deserialization of a FITS file containing a multi-order map.
    ///
    /// The coverage computed contains the cells summing from ``cumul_from`` to ``cumul_to``.
    ///
    /// # Arguments
    ///
    /// * ``cumul_from`` - The cumulative value from which cells are put in the coverage
    /// * ``cumul_to`` - The cumulative value to which cells are put in the coverage
    /// * ``max_depth`` -  the largest depth of the output MOC, which must be larger or equals to the largest
    ///   depth in the `uniq` values
    /// * `asc`: cumulative value computed from lower to highest densities instead of from highest to lowest
    /// * `strict`: (sub-)cells overlapping the `cumul_from` or `cumul_to` values are not added
    /// * `no_split`: cells overlapping the `cumul_from` or `cumul_to` values are not recursively split
    /// * `reverse_decent`: perform the recursive decent from the highest cell number to the lowest (to be compatible with Aladin)
    ///
    /// # Info
    /// 
    /// We expect the FITS file to be a BINTABLE containing a multi-order map.
    /// In this non-flexible approach, we expect the BINTABLE extension to contains:
    /// 
    /// ```bash
    /// XTENSION= 'BINTABLE'           / binary table extension                         
    /// BITPIX  =                    8 / array data type                                
    /// NAXIS   =                    2 / number of array dimensions 
    /// AXIS1  =                    ?? / length of dimension 1                          
    /// NAXIS2  =                   ?? / length of dimension 2                          
    /// PCOUNT  =                    0 / number of group parameters                     
    /// GCOUNT  =                    1 / number of groups                               
    /// TFIELDS =                   xx / number of table fields
    /// TTYPE1  = 'UNIQ    '                                                            
    /// TFORM1  = 'K       '                                                            
    /// TTYPE2  = 'PROBDENSITY'                                                         
    /// TFORM2  = 'D       '                                                            
    /// TUNIT2  = 'sr-1    ' 
    /// ...
    /// MOC     =                    T                                                  
    /// PIXTYPE = 'HEALPIX '           / HEALPIX pixelisation                           
    /// ORDERING= 'NUNIQ   '           / Pixel ordering scheme: RING, NESTED, or NUNIQ  
    /// COORDSYS= 'C       '           / Ecliptic, Galactic or Celestial (equatorial)   
    /// MOCORDER=                   xx / MOC resolution (best order) 
    /// ...
    /// END
    /// ```
    #[pyfn(m, "spatial_moc_from_multiordermap_fits_file")]
    fn spatial_moc_from_multiordermap_fits_file(
        py: Python, 
        path: String,
        cumul_from: f64,
        cumul_to: f64,
        asc: bool,
        strict: bool,
        no_split: bool,
        reverse_decent: bool
    ) -> PyResult<Py<PyArray2<u64>>> {
        use std::fs::File;
        use std::io::BufReader;
        use moc::deser::fits;
        
        let file = File::open(&path)?;
        let reader = BufReader::new(file);
        let ranges = fits::multiordermap::from_fits_multiordermap(
            reader,
            cumul_from, cumul_to,
            asc, strict, no_split, reverse_decent
        ).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))?;
        let result: Array2<u64> = mocranges_to_array2(ranges.into_moc_ranges());
        Ok(result.into_pyarray(py).to_owned())
    }
    
    /// Deserialize a spatial MOC from a FITS file.
    ///
    /// # Arguments
    ///
    /// * ``path`` - The file path
    #[pyfn(m, "spatial_moc_from_fits_file")]
    fn spatial_moc_from_fits_file(py: Python, path: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = spatial_coverage::from_fits_file(path)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a spatial MOC from an ASCII file.
    ///
    /// # Arguments
    ///
    /// * ``path`` - The file path
    #[pyfn(m, "spatial_moc_from_ascii_file")]
    fn spatial_moc_from_ascii_file(py: Python, path: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = spatial_coverage::from_ascii_file(path)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a spatial MOC from a ASCII string.
    ///
    /// # Arguments
    ///
    /// * ``ascii`` - The json string
    #[pyfn(m, "spatial_moc_from_ascii_str")]
    fn spatial_moc_from_ascii_str(py: Python, ascii: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = spatial_coverage::from_ascii_str(ascii)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a spatial MOC from a JSON file.
    ///
    /// # Arguments
    ///
    /// * ``path`` - The file path
    #[pyfn(m, "spatial_moc_from_json_file")]
    fn spatial_moc_from_json_file(py: Python, path: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = spatial_coverage::from_json_file(path)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a spatial MOC from a JSON string.
    ///
    /// # Arguments
    ///
    /// * ``json`` - The json string
    #[pyfn(m, "spatial_moc_from_json_str")]
    fn spatial_moc_from_json_str(py: Python, json: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = spatial_coverage::from_json_str(json)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }




    /// Serialize a time MOC into a FITS file.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    /// * ``path`` - The file path
    #[pyfn(m, "time_moc_to_fits_file")]
    fn time_moc_to_fits_file(
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
        path: String,
    ) -> PyResult<()> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_time_ranges_from_py_uncheked(ranges);
        temporal_coverage::to_fits_file(depth, ranges, path)
          .map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Serialize a time MOC into an ASCII file.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    /// * ``path`` - The file path
    #[pyfn(m, "time_moc_to_ascii_file")]
    fn time_moc_to_ascii_file(
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
        path: String,
    ) -> PyResult<()> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_time_ranges_from_py_uncheked(ranges);
        temporal_coverage::to_ascii_file(depth, ranges, path)
          .map_err(exceptions::PyIOError::new_err)
    }

    /// Serialize a time MOC into a ASCII string.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    #[pyfn(m, "time_moc_to_ascii_str")]
    fn time_moc_to_ascii_str(
        py: Python,
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
    ) -> Py<PyString> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_time_ranges_from_py_uncheked(ranges);
        PyString::new(py,&temporal_coverage::to_ascii_str(depth, ranges)).into()
    }

    /// Serialize a time MOC into a JSON file.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    /// * ``path`` - The file path
    #[pyfn(m, "time_moc_to_json_file")]
    fn time_moc_to_json_file(
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
        path: String,
    ) -> PyResult<()> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_time_ranges_from_py_uncheked(ranges);
        temporal_coverage::to_json_file(depth, ranges, path)
          .map_err(exceptions::PyIOError::new_err)
    }

    /// Serialize a time MOC into a JSON string.
    ///
    /// # Arguments
    ///
    /// * `depth``` - The depth of the MOC (needed to support the case in which there is no cell
    ///               at the deepest level, in which case the computed depth will not be deep enough)
    /// * ``ranges`` - The list of time ranges to serialize.
    #[pyfn(m, "time_moc_to_json_str")]
    fn time_moc_to_json_str(
        py: Python,
        depth: u8,
        ranges: PyReadonlyArray2<u64>,
    ) -> Py<PyString> {
        let ranges = ranges.as_array().to_owned();
        let ranges = coverage::create_time_ranges_from_py_uncheked(ranges);
        PyString::new(py,&temporal_coverage::to_json_str(depth, ranges)).into()
    }

    /// Deserialize a time MOC from a FITS file.
    ///
    /// # Arguments
    ///
    /// * ``path`` - The file path
    #[pyfn(m, "time_moc_from_fits_file")]
    fn time_moc_from_fits_file(py: Python, path: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = temporal_coverage::from_fits_file(path)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a time MOC from an ASCII file.
    ///
    /// # Arguments
    ///
    /// * ``path`` - The file path
    #[pyfn(m, "time_moc_from_ascii_file")]
    fn time_moc_from_ascii_file(py: Python, path: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = temporal_coverage::from_ascii_file(path)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a time MOC from a ASCII string.
    ///
    /// # Arguments
    ///
    /// * ``ascii`` - The json string
    #[pyfn(m, "time_moc_from_ascii_str")]
    fn time_moc_from_ascii_str(py: Python, ascii: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = temporal_coverage::from_ascii_str(ascii)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a time MOC from a JSON file.
    ///
    /// # Arguments
    ///
    /// * ``path`` - The file path
    #[pyfn(m, "time_moc_from_json_file")]
    fn time_moc_from_json_file(py: Python, path: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = temporal_coverage::from_json_file(path)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Deserialize a time MOC from a JSON string.
    ///
    /// # Arguments
    ///
    /// * ``json`` - The json string
    #[pyfn(m, "time_moc_from_json_str")]
    fn time_moc_from_json_str(py: Python, json: String) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = temporal_coverage::from_json_str(json)?;
        let result: Array2<u64> = mocranges_to_array2(ranges);
        Ok(result.into_pyarray(py).to_owned())
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
    /// * ``depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
    #[pyfn(m, "hpx_coverage_degrade")]
    fn hpx_coverage_degrade(
        py: Python,
        ranges: PyReadonlyArray2<u64>,
        depth: u8,
    ) -> PyResult<Py<PyArray2<u64>>> {
        coverage_degrade(py, ranges, depth, coverage::create_hpx_ranges_from_py_unchecked)
    }

    /// Expand the spatial coverage adding an external edge of max_depth pixels.
    ///
    /// # Arguments
    ///
    /// * ``max_depth`` - The MOC depth.
    /// * ``ranges`` - The spatial coverage ranges of max depth to be expanded.
    ///
    /// # Errors
    ///
    /// * ``depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
    #[pyfn(m, "hpx_coverage_expand")]
    fn hpx_coverage_expand(
        py: Python,
        max_depth: u8,
        ranges: PyReadonlyArray2<u64>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = ranges.as_array().to_owned();
        let coverage = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        let result = spatial_coverage::expand(max_depth, coverage);
        let result = mocranges_to_array2(result);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Contract the spatial coverage removing an internal edge of max_depth pixels.
    ///
    /// # Arguments
    ///
    /// * ``max_depth`` - The MOC depth.
    /// * ``ranges`` - The spatial coverage ranges of max depth to be contracted.
    ///
    /// # Errors
    ///
    /// * ``depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
    #[pyfn(m, "hpx_coverage_contract")]
    fn hpx_coverage_contract(
        py: Python,
        max_depth: u8,
        ranges: PyReadonlyArray2<u64>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = ranges.as_array().to_owned();
        let coverage = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        let result = spatial_coverage::contract(max_depth, coverage);
        let result = mocranges_to_array2(result);
        Ok(result.into_pyarray(py).to_owned())
    }

    /// Count the number of disjoint MOC this MOC contains.
    ///
    /// # Arguments
    ///
    /// * ``max_depth`` - The MOC depth.
    /// * ``include_indirect_neighbours`` -
    ///     if `false`, only consider  cells having a common edge as been part of a same MOC
    ///     if `true`, also consider cells having a common vertex as been part of the same MOC
    /// * ``ranges`` - The spatial coverage ranges of max depth to be split.
    ///
    /// # Errors
    ///
    /// * ``depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
    #[pyfn(m, "hpx_coverage_split_count")]
    fn hpx_coverage_split_count(
        max_depth: u8,
        include_indirect_neighbours: bool,
        ranges: PyReadonlyArray2<u64>,
    ) -> u32 {
        let ranges = ranges.as_array().to_owned();
        let coverage = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        let moc = RangeMOC::<u64, Hpx<u64>>::new(max_depth, coverage);
        moc.split_into_joint_mocs(include_indirect_neighbours).len() as u32
    }

    /// Split the input MOC into disjoint MOCs.
    ///
    /// # Arguments
    ///
    /// * ``max_depth`` - The MOC depth.
    /// * ``include_indirect_neighbours`` -
    ///     if `false`, only consider  cells having a common edge as been part of a same MOC
    ///     if `true`, also consider cells having a common vertex as been part of the same MOC
    /// * ``ranges`` - The spatial coverage ranges of max depth to be split.
    ///
    /// # Errors
    ///
    /// * ``depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
    #[pyfn(m, "hpx_coverage_split")]
    fn hpx_coverage_split(
        py: Python,
        max_depth: u8,
        include_indirect_neighbours: bool,
        ranges: PyReadonlyArray2<u64>,
    ) -> PyResult<Py<PyList>> {
        let ranges = ranges.as_array().to_owned();
        let coverage = coverage::create_hpx_ranges_from_py_unchecked(ranges);
        let moc = RangeMOC::<u64, Hpx<u64>>::new(max_depth, coverage);
        let mocs: Vec<Py<PyArray2<u64>>> = moc.split_into_joint_mocs(include_indirect_neighbours)
          .drain(..)
          .map(|cell_moc| 
            vec_range_to_array2(
                cell_moc.into_cell_moc_iter().ranges().collect()
            ).into_pyarray(py).to_owned().into()
          ).collect();
        PyList::new(py, mocs).extract()
    }
    
    
    
    /// Degrade a time coverage to a specific depth.
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The time coverage ranges to degrade.
    /// * ``depth`` - The depth to degrade the time coverage to.
    ///
    /// # Errors
    ///
    /// * ``depth`` is not comprised in `[0, Time::<T>::MAX_DEPTH] = [0, 62]`
    #[pyfn(m, "time_coverage_degrade")]
    fn time_coverage_degrade(
        py: Python,
        ranges: PyReadonlyArray2<u64>,
        depth: u8,
    ) -> PyResult<Py<PyArray2<u64>>> {
        coverage_degrade(py, ranges, depth, coverage::create_time_ranges_from_py_uncheked)
    }


    /// Make a generic coverage consistent
    ///
    /// # Infos
    ///
    /// This is an internal method whose purpose is not to be called
    /// by an user. It is called inside of the `mocpy.IntervalSet` class.
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The coverage ranges to make consistent.
    #[pyfn(m, "coverage_merge_gen_intervals")]
    fn coverage_merge_gen_intervals(
        py: Python,
        ranges: PyReadonlyArray2<u64>
    ) -> PyResult<Py<PyArray2<u64>>> {
        let ranges = ranges.as_array().to_owned();

        let coverage = coverage::build_ranges_from_py(ranges);

        let result: Array2<u64> = ranges_to_array2(coverage);
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
    /// * ``min_depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
    #[pyfn(m, "coverage_merge_hpx_intervals")]
    fn coverage_merge_hpx_intervals(
        py: Python,
        ranges: PyReadonlyArray2<u64>,
        min_depth: i8,
    ) -> PyResult<Py<PyArray2<u64>>> {
        coverage_merge_intervals(py, ranges, min_depth, coverage::build_hpx_ranges_from_py)
    }

    /// Make a time coverage consistent
    ///
    /// # Infos
    ///
    /// This is an internal method whose purpose is not to be called
    /// by an user. It is called inside of the `mocpy.IntervalSet` class.
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The time coverage ranges to make consistent.
    /// * ``min_depth`` - A minimum depth. This argument is optional.
    ///   A min depth means that we do not want any cells to be
    ///   of depth < to ``min_depth``. This argument is used for example for
    ///   plotting a time coverage. All time cells of depth < 2 are splitted
    ///   into cells of depth 2.
    ///
    /// # Errors
    ///
    /// * ``min_depth`` is not comprised in `[0, Time::<T>::MAX_DEPTH] = [0, 62]`
    #[pyfn(m, "coverage_merge_time_intervals")]
    fn coverage_merge_time_intervals(
        py: Python,
        ranges: PyReadonlyArray2<u64>,
        min_depth: i8,
    ) -> PyResult<Py<PyArray2<u64>>> {
        coverage_merge_intervals(py, ranges, min_depth, coverage::build_time_ranges_from_py)
    }

    /// Compute the depth of a spatial coverage
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The input coverage.
    #[pyfn(m, "hpx_coverage_depth")]
    fn hpx_coverage_depth(py: Python, ranges: PyReadonlyArray2<u64>) -> u8 {
        coverage_depth(py, ranges, coverage::create_hpx_ranges_from_py_unchecked)
    }

    /// Compute the depth of a time coverage
    ///
    /// # Arguments
    ///
    /// * ``ranges`` - The input coverage.
    #[pyfn(m, "time_coverage_depth")]
    fn time_coverage_depth(py: Python, ranges: PyReadonlyArray2<u64>) -> u8 {
        coverage_depth(py, ranges, coverage::create_time_ranges_from_py_uncheked)
    }

    fn coverage_depth<Q, F>(_py: Python, ranges: PyReadonlyArray2<u64>, to_moc_ranges: F) -> u8
        where
          Q: MocQty<u64>,
          F: Fn(Array<u64, Ix2>) -> MocRanges<u64, Q>
    {
        let ranges = ranges.as_array().to_owned();
        let coverage = to_moc_ranges(ranges);
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
            Array::zeros((0, 2))
        } else {
            let shape = (ranges.shape()[0], 1);
            let start = ranges.into_shape(shape).unwrap();
            let mut ranges: Vec<Range<u64>> = Vec::with_capacity(start.len());
            for uniq in start {
                ranges.push(uniq..uniq + 1)
            }
            ranges.sort_by(|a, b| a.start.cmp(&b.start));
            let nested_coverage = spatial_coverage::to_nested(HpxUniqRanges::new_unchecked(ranges));
            mocranges_to_array2(nested_coverage)
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
        use moc::moc::range::RangeMOC;
        use moc::moc::{RangeMOCIterator, RangeMOCIntoIterator};
        
        let ranges = ranges.as_array().to_owned();

        let result: Array1<u64> = if ranges.is_empty() {
            Array::zeros((0,))
        } else {
            let nested_coverage = coverage::create_hpx_ranges_from_py_unchecked(ranges);

            // let uniq_coverage = nested_coverage.into_hpx_uniq();
            // uniq_ranges_to_array1(uniq_coverage)
            
            let mut v: Vec<u64> = RangeMOC::new(29, nested_coverage).into_range_moc_iter()
              .cells()
              .map(|cell| cell.uniq_hpx())
              .collect();
            v.sort_unstable();
            v.into()
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
    /// # WARNING
    /// * using `f64`, it is not precise to the microsecond,
    ///   use `from_time_ranges_in_microsec_since_jd_origin` instead.
    ///
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

        let coverage: Array2<u64> = temporal_coverage::from_time_ranges(min_times, max_times)?;

        Ok(coverage.into_pyarray(py).to_owned())
    }

    /// Create a temporal coverage from a list of time ranges expressed in microseconds since
    /// jd origin.
    ///
    /// # Arguments
    ///
    /// * ``min_times`` - The list of inf bounds of the time ranges expressed in microseconds since
    ///    jd origin.
    /// * ``max_times`` - The list of sup bounds of the time ranges expressed in microseconds since
    ///    jd origin.
    ///
    /// # Errors
    ///
    /// * If the number of ``min_times`` and ``max_times`` do not match.
    #[pyfn(m, "from_time_ranges_in_microsec_since_jd_origin")]
    fn from_time_ranges_in_microsec_since_jd_origin(
        py: Python,
        min_times: PyReadonlyArray1<u64>,
        max_times: PyReadonlyArray1<u64>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let min_times = min_times.as_array().to_owned();
        let max_times = max_times.as_array().to_owned();

        let coverage: Array2<u64> = temporal_coverage::from_time_ranges_in_microsec_since_jd_origin(min_times, max_times)?;

        Ok(coverage.into_pyarray(py).to_owned())
    }

    /// Flatten HEALPix cells to a specific depth
    ///
    /// # Arguments
    ///
    /// * ``data`` - The spatial coverage
    /// * ``depth`` - The depth to flatten the coverage to.
    #[pyfn(m, "flatten_pixels")]
    fn flatten_hpx_pixels(py: Python, data: PyReadonlyArray2<u64>, depth: u8) -> Py<PyArray1<u64>> {
        let data = data.as_array().to_owned();

        let result = coverage::flatten_hpx_pixels(data, depth);

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
        depth: PyReadonlyArray1<u8>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let pixels = pixels.as_array().to_owned();
        let depth = depth.as_array().to_owned();

        let result = spatial_coverage::from_healpix_cells(pixels, depth)?;

        Ok(result.into_pyarray(py).to_owned())
    }


    /// Create a spatial coverage from an HEALPix map, i.e. from a list of HEALPix cell indices
    /// at the same depth.
    ///
    /// # Arguments
    ///
    /// * ``pixels`` - A set of HEALPix cell indices
    /// * ``depth`` - The depths of each HEALPix cell indices
    ///
    /// # Precondition
    ///
    /// * ``depth`` is a value in the range `[0, <T>::MAXDEPTH] = [0, 29]`
    /// * ``pixels`` contains values in the range `[0, 12*4**(depth)]`
    #[pyfn(m, "from_healpix_cells")]
    fn from_healpix_map(
        py: Python,
        pixels: PyReadonlyArray1<u64>,
        depth: PyReadonlyArray1<u8>,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let pixels = pixels.as_array().to_owned();
        let depth = depth.as_array().to_owned();

        let result = spatial_coverage::from_healpix_cells(pixels, depth)?;

        Ok(result.into_pyarray(py).to_owned())
    }

    /// Create a spatial coverage from a given ring.
    ///
    /// # Arguments
    ///
    /// * ``ra_deg`` - longitude of the center of the ring, in degrees
    /// * ``dec_deg`` - latitude of the center of the ring, in degrees
    /// * ``r_int_deg`` - Internal radius of the ring, in degrees
    /// * ``r_ext_deg`` - External radius of the ring, in degrees
    /// * ``depth`` - The depths of the expected MOC
    /// * ``delta_depth`` - parameter controlling the approximation (typical value: 2)
    ///
    /// # Errors
    ///
    /// If one of the following conditions is not met:
    ///
    /// * ``depth`` contains values in the range `[0, <T>::MAXDEPTH] = [0, 29]`
    /// * ``r_int_deg`` contains values in the range `[0, 180]`
    /// * ``r_ext_deg`` contains values in the range `[0, 180]`
    /// * ``r_ext_deg > r_int_deg``
    ///
    #[pyfn(m, "from_ring")]
    fn from_ring(
        py: Python,
        lon_deg: f64,
        lat_deg: f64,
        r_int_deg: f64,
        r_ext_deg: f64,
        depth: u8,
        delta_depth: u8,
    ) -> PyResult<Py<PyArray2<u64>>> {
        let moc_ranges = RangeMOC::from_ring(
            lon_deg.to_radians(),
            lat_deg.to_radians(),
            r_int_deg.to_radians(),
            r_ext_deg.to_radians(),
            depth,
            delta_depth
        ).into_moc_ranges();
        let result = mocranges_to_array2(moc_ranges);
        Ok(result.into_pyarray(py).to_owned())
    }

   Ok(())
}
