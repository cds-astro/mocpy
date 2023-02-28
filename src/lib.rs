
use std::ops::Range;

use ndarray::Array;
use numpy::{
  Ix2, Ix3,
  IntoPyArray,
  PyArray1, PyArray2, PyArray3, PyArrayDyn,
  PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArrayDyn
};
use pyo3::{
  types::{PyBytes, PyTuple},
  prelude::{pymodule, Py, PyModule, PyResult, Python},
  exceptions::{PyIOError, PyValueError},
};

use moc::{
  utils,
  storage::u64idx::U64MocStore
};

#[pymodule]
fn mocpy(_py: Python, m: &PyModule) -> PyResult<()> {
  #[pyfn(m)]
  fn usize_n_bits() -> u8 {
    (std::mem::size_of::<usize>() as u8) << 3
  }

  /// Make a new spatial coverage from ranges (at max order) and a depth
  ///
  ///
  /// # Arguments
  ///
  /// * ``depth`` - The depth of the coverage between `[0, <u64>::MAXDEPTH] = [0, 29]`
  /// * ``ranges`` - The coverage ranges.
  #[pyfn(m)]
  fn from_hpx_ranges(
    _py: Python,
    depth: u8,
    ranges: PyReadonlyArray2<u64>,
  ) -> PyResult<usize> {
    let it = ranges.as_array().into_iter().cloned();
    struct RangeIt<T: Iterator<Item=u64>> {
      it: T
    }
    impl<T: Iterator<Item=u64>> Iterator for RangeIt<T> {
      type Item = Range<u64>;
      fn next(&mut self) -> Option<Self::Item> {
        if let (Some(start), Some(end)) = (self.it.next(), self.it.next()) {
          Some(start..end)
        } else {
          None
        }
      }
    }
    U64MocStore::get_global_store()
      .from_hpx_ranges(depth, RangeIt { it }, None)
      .map_err(PyValueError::new_err)
  }

  /// Create a temporal coverage from a list of time ranges expressed in microseconds since
  /// jd origin.
  ///
  /// # Arguments
  ///
  /// * ``depth`` - depth of the MOC
  /// * ``ranges``: PyReadonlyArray2<u64>,
  #[pyfn(m)]
  fn from_time_ranges_array2(
    depth: u8,
    ranges: PyReadonlyArray2<u64>,
  ) -> PyResult<usize> {
    let it = ranges.as_array().into_iter().cloned();
    struct RangeIt<T: Iterator<Item=u64>> {
      it: T
    }
    impl<T: Iterator<Item=u64>> Iterator for RangeIt<T> {
      type Item = Range<u64>;
      fn next(&mut self) -> Option<Self::Item> {
        if let (Some(start), Some(end)) = (self.it.next(), self.it.next()) {
          Some(start..end)
        } else {
          None
        }
      }
    }
    U64MocStore::get_global_store()
      .from_microsec_ranges_since_jd0(depth, RangeIt { it })
      .map_err(PyValueError::new_err)
  }


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
  #[pyfn(m)]
  fn from_lonlat(
    depth: u8,
    lon_deg: PyReadonlyArrayDyn<f64>,
    lat_deg: PyReadonlyArrayDyn<f64>,
  ) -> PyResult<usize> {
    let lon = lon_deg.as_array().into_iter().cloned();
    let lat = lat_deg.as_array().into_iter().cloned();
    U64MocStore::get_global_store()
      .from_coo(depth, lon.zip(lat))
      .map_err(PyValueError::new_err)
  }

  #[pyfn(m)]
  fn from_cone(
    lon_deg: f64,
    lat_deg: f64,
    radius_deg: f64,
    depth: u8,
    delta_depth: u8
  ) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .from_cone(lon_deg, lat_deg, radius_deg, depth, delta_depth)
      .map_err(PyValueError::new_err)
  }

  /// Create a spatial coverage from a given ring.
  ///
  /// # Arguments
  ///
  /// * ``lon_deg`` - longitude of the center of the ring, in degrees
  /// * ``lat_deg`` - latitude of the center of the ring, in degrees
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
  #[pyfn(m)]
  fn from_ring(
    lon_deg: f64,
    lat_deg: f64,
    r_int_deg: f64,
    r_ext_deg: f64,
    depth: u8,
    delta_depth: u8,
  ) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .from_ring(
        lon_deg,
        lat_deg,
        r_int_deg,
        r_ext_deg,
        depth,
        delta_depth,
      )
      .map_err(PyValueError::new_err)
  }

  #[pyfn(m)]
  pub fn from_elliptical_cone(
    lon_deg: f64,
    lat_deg: f64,
    a_deg: f64,
    b_deg: f64,
    pa_deg: f64,
    depth: u8,
    delta_depth: u8,
  ) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .from_elliptical_cone(lon_deg, lat_deg, a_deg, b_deg, pa_deg, depth, delta_depth)
      .map_err(PyValueError::new_err)
  }

  #[pyfn(m)]
  pub fn from_polygon(
    lon_deg: PyReadonlyArrayDyn<f64>,
    lat_deg: PyReadonlyArrayDyn<f64>,
    complement: bool,
    depth: u8,
  ) -> PyResult<usize> {
    let lon = lon_deg.as_array().into_iter().cloned();
    let lat = lat_deg.as_array().into_iter().cloned();
    U64MocStore::get_global_store()
      .from_polygon(lon.zip(lat), complement, depth)
      .map_err(PyValueError::new_err)
  }

  /// Create a 1D spatial coverage from a list of uniq cells each associated with a value.
  ///
  /// The coverage computed contains the cells summing from ``cumul_from`` to ``cumul_to``.
  ///
  /// # Arguments
  ///
  /// * ``max_depth`` -  the largest depth of the output MOC, which must be larger or equals to the largest
  ///                    depth in the `uniq` values
  /// * ``uniq`` - Uniq HEALPix indices
  /// * ``values`` - Array containing the values associated to each cell.
  /// * ``values_are_densities`` - sum of all values equals 1, the values are densities (doe not depend on cell size)
  /// * ``cumul_from`` - The cumulative value from which cells are put in the coverage
  /// * ``cumul_to`` - The cumulative value to which cells are put in the coverage
  /// * `asc`: cumulative value computed from lower to highest densities instead of from highest to lowest
  /// * `strict`: (sub-)cells overlapping the `cumul_from` or `cumul_to` values are not added
  /// * `no_split`: cells overlapping the `cumul_from` or `cumul_to` values are not recursively split
  /// * `reverse_decent`: perform the recursive decent from the highest cell number to the lowest (to be compatible with Aladin)
  ///
  #[pyfn(m)]
  fn from_valued_hpx_cells(
    max_depth: u8,
    uniq: PyReadonlyArrayDyn<u64>,
    values: PyReadonlyArrayDyn<f64>,
    values_are_densities: bool,
    cumul_from: f64,
    cumul_to: f64,
    asc: bool,
    strict: bool,
    no_split: bool,
    reverse_decent: bool,
  ) -> PyResult<usize> {
    let uniq = uniq.as_array();
    let values = values.as_array();
    if uniq.len() != values.len() {
      return Err(PyValueError::new_err("`uniq` and values do not have the same size."));
    }
    let uniq_vals = uniq.into_iter().cloned().zip(values.into_iter().cloned());
    U64MocStore::get_global_store()
      .from_valued_cells(
        max_depth,
        values_are_densities,
        cumul_from,
        cumul_to,
        asc,
        !strict,
        !no_split,
        reverse_decent,
        uniq_vals,
      )
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn from_time_lonlat_approx(
    times: PyReadonlyArrayDyn<f64>,
    d1: u8,
    lon: PyReadonlyArrayDyn<f64>,
    lat: PyReadonlyArrayDyn<f64>,
    d2: u8,
  ) -> PyResult<usize> {
    let times = times.as_array().to_owned().into_raw_vec();
    let lon = lon.as_array().to_owned().into_raw_vec();
    let lat = lat.as_array().to_owned().into_raw_vec();

    U64MocStore::get_global_store()
      .create_from_times_positions_approx(times, lon, lat, d1, d2)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn from_time_lonlat(
    times: PyReadonlyArrayDyn<u64>,
    d1: u8,
    lon: PyReadonlyArrayDyn<f64>,
    lat: PyReadonlyArrayDyn<f64>,
    d2: u8,
  ) -> PyResult<usize> {
    let times = times.as_array().to_owned().into_raw_vec();
    let lon = lon.as_array().to_owned().into_raw_vec();
    let lat = lat.as_array().to_owned().into_raw_vec();

    U64MocStore::get_global_store()
      .create_from_times_positions(times, lon, lat, d1, d2)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn from_time_ranges_lonlat_approx(
    times_min: PyReadonlyArrayDyn<f64>,
    times_max: PyReadonlyArrayDyn<f64>,
    d1: u8,
    lon: PyReadonlyArrayDyn<f64>,
    lat: PyReadonlyArrayDyn<f64>,
    d2: u8,
  ) -> PyResult<usize> {
    let times_min = times_min.as_array().to_owned().into_raw_vec();
    let times_max = times_max.as_array().to_owned().into_raw_vec();
    let lon = lon.as_array().to_owned().into_raw_vec();
    let lat = lat.as_array().to_owned().into_raw_vec();

    U64MocStore::get_global_store()
      .create_from_time_ranges_positions_approx(times_min, times_max, d1, lon, lat, d2)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn from_time_ranges_lonlat(
    times_min: PyReadonlyArrayDyn<u64>,
    times_max: PyReadonlyArrayDyn<u64>,
    d1: u8,
    lon: PyReadonlyArrayDyn<f64>,
    lat: PyReadonlyArrayDyn<f64>,
    d2: u8,
  ) -> PyResult<usize> {
    let times_min = times_min.as_array().to_owned().into_raw_vec();
    let times_max = times_max.as_array().to_owned().into_raw_vec();
    let lon = lon.as_array().to_owned().into_raw_vec();
    let lat = lat.as_array().to_owned().into_raw_vec();

    U64MocStore::get_global_store()
      .create_from_time_ranges_positions(times_min, times_max, d1, lon, lat, d2)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn from_time_ranges_spatial_coverages_approx(
    times_min: PyReadonlyArrayDyn<f64>,
    times_max: PyReadonlyArrayDyn<f64>,
    d1: u8,
    spatial_coverages: PyReadonlyArrayDyn<usize>,
  ) -> PyResult<usize> {
    let times_min = times_min.as_array().to_owned().into_raw_vec();
    let times_max = times_max.as_array().to_owned().into_raw_vec();
    if times_min.len() != times_max.len() {
      return Err(PyValueError::new_err("`times_min` and `times_max` do not have the same size."));
    }
    let spatial_coverage_indices = spatial_coverages.as_array().to_owned().into_raw_vec();
    if times_min.len() != spatial_coverage_indices.len() {
      return Err(PyValueError::new_err("`times` and `spatial indices` do not have the same size."));
    }
    U64MocStore::get_global_store()
      .from_time_ranges_spatial_coverages_in_store_approx(times_min, times_max, d1, spatial_coverage_indices)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn from_time_ranges_spatial_coverages(
    times_min: PyReadonlyArrayDyn<u64>,
    times_max: PyReadonlyArrayDyn<u64>,
    d1: u8,
    spatial_coverages: PyReadonlyArrayDyn<usize>,
  ) -> PyResult<usize> {
    let times_min = times_min.as_array().to_owned().into_raw_vec();
    let times_max = times_max.as_array().to_owned().into_raw_vec();
    if times_min.len() != times_max.len() {
      return Err(PyValueError::new_err("`times_min` and `times_max` do not have the same size."));
    }
    let spatial_coverage_indices = spatial_coverages.as_array().to_owned().into_raw_vec();
    if times_min.len() != spatial_coverage_indices.len() {
      return Err(PyValueError::new_err("`times` and `spatial indices` do not have the same size."));
    }

    U64MocStore::get_global_store()
      .from_time_ranges_spatial_coverages_in_store(times_min, times_max, d1, spatial_coverage_indices)
      .map_err(PyValueError::new_err)
  }


  #[pyfn(m)]
  fn project_on_first_dim(smoc_index: usize, stmoc_index: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .space_fold(smoc_index, stmoc_index)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn project_on_second_dim(tmoc_index: usize, stmoc_index: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .time_fold(tmoc_index, stmoc_index)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_from_fits_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_stmoc_from_fits_file(path)
      .map_err(PyIOError::new_err)
  }

  #[pyfn(m)]
  fn stmoc_from_fits_raw_bytes(raw_bytes: &[u8]) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_stmoc_from_fits_buff(raw_bytes)
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_from_ascii_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_stmoc_from_ascii_file(path)
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_from_json_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_stmoc_from_json_file(path)
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_from_ascii_str(_py: Python, ascii: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_stmoc_from_ascii(&ascii)
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_from_json_str(_py: Python, json: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_stmoc_from_json(&json)
      .map_err(PyIOError::new_err)
  }


  /// Drop the MOC at the given index.
  ///
  /// This method is automatically called by the Python garbage collector.
  #[pyfn(m)]
  fn copy(index: usize) -> PyResult<()> {
    U64MocStore::get_global_store()
      .copy(index)
      .map_err(PyIOError::new_err)
  }

  /// Drop the MOC at the given index.
  ///
  /// This method is automatically called by the Python garbage collector.
  #[pyfn(m)]
  fn drop(index: usize) -> PyResult<()> {
    U64MocStore::get_global_store()
      .drop(index)
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_depth(_py: Python, index: usize) -> PyResult<(u8, u8)> {
    U64MocStore::get_global_store()
      .get_stmoc_depths(index)
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_min_time_approx(_py: Python, index: usize) -> PyResult<f64> {
    U64MocStore::get_global_store()
      .get_1st_axis_min(index)
      .and_then(|opt| opt.ok_or_else(|| String::from("Empty ST-MOC")))
      .map(|t_min_microsec_since_jd0| (t_min_microsec_since_jd0 as f64) / 86400000000_f64)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_min_time(_py: Python, index: usize) -> PyResult<u64> {
    U64MocStore::get_global_store()
      .get_1st_axis_min(index)
      .and_then(|opt| opt.ok_or_else(|| String::from("Empty ST-MOC")))
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_max_time_approx(_py: Python, index: usize) -> PyResult<f64> {
    U64MocStore::get_global_store()
      .get_1st_axis_max(index)
      .and_then(|opt| opt.ok_or_else(|| String::from("Empty ST-MOC")))
      .map(|t_min_microsec_since_jd0| (t_min_microsec_since_jd0 as f64) / 86400000000_f64)
      .map_err(PyValueError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_max_time(_py: Python, index: usize) -> PyResult<u64> {
    U64MocStore::get_global_store()
      .get_1st_axis_max(index)
      .and_then(|opt| opt.ok_or_else(|| String::from("Empty ST-MOC")))
      .map_err(PyValueError::new_err)
  }


  /// Computes the complement of the given coverage
  ///
  /// # Arguments
  ///
  /// * ``id`` - index of the coverage
  #[pyfn(m)]
  fn complement(_py: Python, id: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .not(id)
      .map_err(PyIOError::new_err)
  }

  /// Perform the union between two coverages of same type.
  ///
  /// # Arguments
  ///
  /// * ``id_left`` - index of the coverage being in the left of the operation.
  /// * ``id_right`` - index of the coverage being in the right of the operation.
  #[pyfn(m)]
  fn union(_py: Python, id_left: usize, id_right: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .or(id_left, id_right)
      .map_err(PyIOError::new_err)
  }

  /// Perform the union between multiple coverages of same type.
  ///
  /// # Arguments
  ///
  /// * ``ids`` - indices of the coverages.
  #[pyfn(m)]
  fn multi_union(_py: Python, ids: PyReadonlyArray1<usize>) -> PyResult<usize> {
    let indices = ids.as_slice()?;
    U64MocStore::get_global_store()
      .multi_union(indices)
      .map_err(PyIOError::new_err)
  }

  /// Perform the intersection between two coverages of same type.
  ///
  /// # Arguments
  ///
  /// * ``id_left`` - index of the coverage being in the left of the operation.
  /// * ``id_right`` - index of the coverage being in the right of the operation.
  #[pyfn(m)]
  fn intersection(_py: Python, id_left: usize, id_right: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .and(id_left, id_right)
      .map_err(PyIOError::new_err)
  }

  /// Perform the intersection between multiple coverages of same type.
  ///
  /// # Arguments
  ///
  /// * ``ids`` - indices of the coverages.
  #[pyfn(m)]
  fn multi_intersection(_py: Python, ids: PyReadonlyArray1<usize>) -> PyResult<usize> {
    let indices = ids.as_slice()?;
    U64MocStore::get_global_store()
      .multi_intersection(indices)
      .map_err(PyIOError::new_err)
  }

  /// Perform the difference between two coverages of same type.
  ///
  /// # Arguments
  ///
  /// * ``id_left`` - index of the coverage being in the left of the operation.
  /// * ``id_right`` - index of the coverage being in the right of the operation.
  #[pyfn(m)]
  fn symmetric_difference(_py: Python, id_left: usize, id_right: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .symmetric_difference(id_left, id_right)
      .map_err(PyIOError::new_err)
  }

  /// Perform the difference between multiple coverages of same type.
  ///
  /// # Arguments
  ///
  /// * ``ids`` - indices of the coverages.
  #[pyfn(m)]
  fn multi_symmetric_difference(_py: Python, ids: PyReadonlyArray1<usize>) -> PyResult<usize> {
    let indices = ids.as_slice()?;
    U64MocStore::get_global_store()
      .multi_symmetric_difference(indices)
      .map_err(PyIOError::new_err)
  }

  /// Perform the difference between two coverages of same type.
  ///
  /// # Arguments
  ///
  /// * ``id_left`` - index of the coverage being in the left of the operation.
  /// * ``id_right`` - index of the coverage being in the right of the operation.
  #[pyfn(m)]
  fn difference(_py: Python, id_left: usize, id_right: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .difference(id_left, id_right)
      .map_err(PyIOError::new_err)
  }

  /// Check the equality between two coverages
  ///
  /// # Arguments
  ///
  /// * ``id_left`` - index of the coverage being in the left of the operation.
  /// * ``id_right`` - index of the coverage being in the right of the operation.
  #[pyfn(m)]
  fn check_eq(_py: Python, id_left: usize, id_right: usize) -> PyResult<bool> {
    U64MocStore::get_global_store()
      .eq(id_left, id_right)
      .map_err(PyIOError::new_err)
  }

  /// Checks whether a coverage is empty.
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage to check the emptiness.
  #[pyfn(m)]
  fn is_empty(_py: Python, index: usize) -> PyResult<bool> {
    U64MocStore::get_global_store()
      .is_empty(index)
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn coverage_2d_contains_approx(
    py: Python,
    index: usize,
    times: PyReadonlyArrayDyn<f64>,
    lon: PyReadonlyArrayDyn<f64>,
    lat: PyReadonlyArrayDyn<f64>
  ) -> PyResult<Py<PyArrayDyn<bool>>> {
    let time_shape = times.shape().to_vec();
    let lon_shape = lon.shape();
    let lat_shape = lat.shape();
    if time_shape != lat_shape {
      return Err(PyValueError::new_err(format!("Time shape different from lon shape: {:?} != {:?}", time_shape, lon_shape)));
    }
    if lon_shape != lat_shape {
      return Err(PyValueError::new_err(format!("Lon shape different from lat shape: {:?} != {:?}", lon_shape, lat_shape)));
    }
    let it_time = times.as_array().into_iter().cloned();
    let it_lon = lon.as_array().into_iter().cloned();
    let it_lat = lat.as_array().into_iter().cloned();
    let it = it_time.zip(it_lon.zip(it_lat));
    U64MocStore::get_global_store()
      .filter_timepos_approx(index, it, |b| b) // in numpy, the mask is reversed (true means do not select)
      .map_err(PyIOError::new_err)
      .and_then(|vec_bool| Array::from_shape_vec(time_shape, vec_bool)
        .map(|a| a.into_pyarray(py).to_owned())
        .map_err(|e| PyValueError::new_err(e.to_string()))
      )
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
  #[pyfn(m)]
  fn coverage_2d_contains(
    py: Python,
    index: usize,
    times: PyReadonlyArrayDyn<u64>,
    lon: PyReadonlyArrayDyn<f64>,
    lat: PyReadonlyArrayDyn<f64>
  ) -> PyResult<Py<PyArray1<bool>>> {
    let it_time = times.as_array().into_iter().cloned();
    let it_lon = lon.as_array().into_iter().cloned();
    let it_lat = lat.as_array().into_iter().cloned();
    // TODO: check all shapes and return an array of same shape!! (see other methods with reshape!)
    let it = it_time.zip(it_lon.zip(it_lat));
    U64MocStore::get_global_store()
      .filter_timepos(index, it, |b| b) // in numpy, the mask is reversed (true means do not select)
      .map(|vec_bool| PyArray1::<bool>::from_vec(py, vec_bool).to_owned())
      .map_err(PyIOError::new_err)
  }

  /// Checks if lon lat coordinates are contained into a Space coverage
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the Space coverage.
  /// * ``lon`` - The longitudes, in degrees.
  /// * ``lat`` - The latitudes, in degrees.
  ///
  /// # Errors
  ///
  /// * If `lon` and `lat` do not have the same length
  #[pyfn(m)]
  fn filter_pos(
    py: Python,
    index: usize,
    lon: PyReadonlyArrayDyn<f64>,
    lat: PyReadonlyArrayDyn<f64>
  ) -> PyResult<Py<PyArrayDyn<bool>>> {
    let lon_shape = lon.shape().to_vec();
    let lat_shape = lat.shape();
    if lon_shape != lat_shape {
      return Err(PyValueError::new_err(format!("Lon shape different from lat shape: {:?} != {:?}", lon_shape, lat_shape)));
    }
    let it_lon = lon.as_array().into_iter().cloned();
    let it_lat = lat.as_array().into_iter().cloned();
    U64MocStore::get_global_store()
      .filter_pos(index, it_lon.zip(it_lat), |b| b)
      .map_err(PyIOError::new_err)
      .and_then(|vec_bool| Array::from_shape_vec(lon_shape, vec_bool)
        .map(|a| a.into_pyarray(py).to_owned())
        .map_err(|e| PyValueError::new_err(e.to_string()))
      )
  }

  /// Checks if a given times are contained into a Time coverage
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the Space coverage.
  /// * ``time`` - The time, in JD.
  ///
  #[pyfn(m)]
  fn filter_time_approx(
    py: Python,
    index: usize,
    times: PyReadonlyArrayDyn<f64>,
  ) -> PyResult<Py<PyArrayDyn<bool>>> {
    let time_shape = times.shape().to_vec();
    let it_time = times.as_array().into_iter().cloned();
    U64MocStore::get_global_store()
      .filter_time_approx(index, it_time, |b| b)
      .map_err(PyIOError::new_err)
      .and_then(|vec_bool| Array::from_shape_vec(time_shape, vec_bool)
        .map(|a| a.into_pyarray(py).to_owned())
        .map_err(|e| PyValueError::new_err(e.to_string()))
      )
  }

  /// Checks if a given times are contained into a Time coverage
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the Space coverage.
  /// * ``time`` - The time, in microsec since JD=0.
  ///
  #[pyfn(m)]
  fn filter_time(
    py: Python,
    index: usize,
    times: PyReadonlyArrayDyn<u64>,
  ) -> PyResult<Py<PyArrayDyn<bool>>> {
    let time_shape = times.shape().to_vec();
    let it_time = times.as_array().into_iter().cloned();
    U64MocStore::get_global_store()
      .filter_time(index, it_time, |b| b)
      .map_err(PyIOError::new_err)
      .and_then(|vec_bool| Array::from_shape_vec(time_shape, vec_bool)
        .map(|a| a.into_pyarray(py).to_owned())
        .map_err(|e| PyValueError::new_err(e.to_string()))
      )
  }

  /// Serialize a coverage into FITS blob.
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage to serialize.
  /// * ``path`` - the path of the output file
  #[pyfn(m)]
  fn to_fits_raw(py: Python, index: usize, pre_v2: bool) -> PyResult<Py<PyBytes>> {
    U64MocStore::get_global_store()
      .to_fits_buff(index, Some(pre_v2))
      .map(move |b| PyBytes::new(py, &b).into())
      .map_err(PyIOError::new_err)
  }

  /// Serialize a coverage into a FITS file
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage to serialize.
  /// * ``path`` - the path of the output file
  #[pyfn(m)]
  fn to_fits_file(index: usize, path: String, pre_v2: bool) -> PyResult<()> {
    U64MocStore::get_global_store()
      .to_fits_file(index, path, Some(pre_v2))
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into an ASCII file.
  ///
  /// # Arguments
  ///
  ///* ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn to_ascii_file(index: usize, path: String) -> PyResult<()> {
    U64MocStore::get_global_store()
      .to_ascii_file(index, path, None)
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into an ASCII file.
  ///
  /// # Arguments
  ///
  ///* ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  /// * ``fold`` - value of the fold parameter (to limit line width)
  #[pyfn(m)]
  fn to_ascii_file_with_fold(index: usize, path: String, fold: usize) -> PyResult<()> {
    U64MocStore::get_global_store()
      .to_ascii_file(index, path, Some(fold))
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into a ASCII string.
  ///
  /// # Arguments
  ///
  /// * ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn to_ascii_str(index: usize) -> PyResult<String> {
    U64MocStore::get_global_store()
      .to_ascii_str(index, None)
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into a ASCII string.
  ///
  /// # Arguments
  ///
  /// * ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  /// * ``fold`` - value of the fold parameter (to limit line width)
  #[pyfn(m)]
  fn to_ascii_str_with_fold(index: usize, fold: usize) -> PyResult<String> {
    U64MocStore::get_global_store()
      .to_ascii_str(index, Some(fold))
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into an JSON file.
  ///
  /// # Arguments
  ///
  /// * ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn to_json_file(index: usize, path: String) -> PyResult<()> {
    U64MocStore::get_global_store()
      .to_json_file(index, path, None)
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into an JSON file.
  ///
  /// # Arguments
  ///
  ///* ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  /// * ``fold`` - value of the fold parameter (to limit line width)
  #[pyfn(m)]
  fn to_json_file_with_fold(index: usize, path: String, fold: usize) -> PyResult<()> {
    U64MocStore::get_global_store()
      .to_json_file(index, path, Some(fold))
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into a JSON string.
  ///
  /// # Arguments
  ///
  /// * ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn to_json_str(index: usize) -> PyResult<String> {
    U64MocStore::get_global_store()
      .to_json_str(index, None)
      .map_err(PyIOError::new_err)
  }

  /// Serialize a MOC into a JSON string.
  ///
  /// # Arguments
  ///
  /// * ``index`` - index in the store of the MOC to be serialized
  /// * ``path`` - The file path
  /// * ``fold`` - value of the fold parameter (to limit line width)
  #[pyfn(m)]
  fn to_json_str_with_fold(index: usize, fold: usize) -> PyResult<String> {
    U64MocStore::get_global_store()
      .to_json_str(index, Some(fold))
      .map_err(PyIOError::new_err)
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
  #[pyfn(m)]
  fn spatial_moc_from_multiordermap_fits_file(
    path: String,
    cumul_from: f64,
    cumul_to: f64,
    asc: bool,
    strict: bool,
    no_split: bool,
    reverse_decent: bool,
  ) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .from_multiordermap_fits_file(
        path, cumul_from, cumul_to, asc, !strict, !no_split, reverse_decent
      )
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a spatial MOC from a FITS file.
  ///
  /// # Arguments
  ///
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn spatial_moc_from_fits_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_smoc_from_fits_file(path)
      .map_err(PyIOError::new_err)
  }

  #[pyfn(m)]
  fn spatial_moc_from_fits_raw_bytes(raw_bytes: &[u8]) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_from_fits_buff(raw_bytes)
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a spatial MOC from an ASCII file.
  ///
  /// # Arguments
  ///
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn spatial_moc_from_ascii_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_smoc_from_ascii_file(path)
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a spatial MOC from a ASCII string.
  ///
  /// # Arguments
  ///
  /// * ``ascii`` - The json string
  #[pyfn(m)]
  fn spatial_moc_from_ascii_str(ascii: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_smoc_from_ascii(ascii.as_str())
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a spatial MOC from a JSON file.
  ///
  /// # Arguments
  ///
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn spatial_moc_from_json_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_smoc_from_json_file(path)
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a spatial MOC from a JSON string.
  ///
  /// # Arguments
  ///
  /// * ``json`` - The json string
  #[pyfn(m)]
  fn spatial_moc_from_json_str(json: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_smoc_from_json(json.as_str())
      .map_err(PyIOError::new_err)
  }



  /// Deserialize a time MOC from a FITS file.
  ///
  /// # Arguments
  ///
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn time_moc_from_fits_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_tmoc_from_fits_file(path)
      .map_err(PyIOError::new_err)
  }

  #[pyfn(m)]
  fn time_moc_from_fits_raw_bytes(raw_bytes: &[u8]) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_tmoc_from_fits_buff(raw_bytes)
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a time MOC from an ASCII file.
  ///
  /// # Arguments
  ///
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn time_moc_from_ascii_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_tmoc_from_ascii_file(path)
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a time MOC from a ASCII string.
  ///
  /// # Arguments
  ///
  /// * ``ascii`` - The json string
  #[pyfn(m)]
  fn time_moc_from_ascii_str(ascii: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_tmoc_from_ascii(ascii.as_str())
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a time MOC from a JSON file.
  ///
  /// # Arguments
  ///
  /// * ``path`` - The file path
  #[pyfn(m)]
  fn time_moc_from_json_file(path: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_tmoc_from_json_file(path)
      .map_err(PyIOError::new_err)
  }

  /// Deserialize a time MOC from a JSON string.
  ///
  /// # Arguments
  ///
  /// * ``json`` - The json string
  #[pyfn(m)]
  fn time_moc_from_json_str(json: String) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .load_tmoc_from_json(json.as_str())
      .map_err(PyIOError::new_err)
  }

  /// Expand the spatial coverage adding an external edge of max_depth pixels
  /// and return the index of the newly created moc.
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage in the store.
  ///
  /// # Errors
  ///
  /// * ``depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
  #[pyfn(m)]
  fn extend(index: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .extend(index)
      .map_err(PyValueError::new_err)
  }

  /// Contract the spatial coverage removing an internal edge of max_depth pixels
  /// and return the index of the newly created moc.
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage in the store.
  ///
  #[pyfn(m)]
  fn contract(index: usize) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .contract(index)
      .map_err(PyValueError::new_err)
  }

  /// Count the number of disjoint MOC this MOC contains.
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage in the store.
  /// * ``include_indirect_neighbours`` -
  ///     if `false`, only consider  cells having a common edge as been part of a same MOC
  ///     if `true`, also consider cells having a common vertex as been part of the same MOC
  ///
  #[pyfn(m)]
  fn split_count(_py: Python, index: usize, include_indirect_neighbours: bool) -> PyResult<u32> {
    if include_indirect_neighbours {
      U64MocStore::get_global_store().split_indirect_count(index)
    } else {
      U64MocStore::get_global_store().split_count(index)
    }.map_err(PyValueError::new_err)
  }

  /// Split the input MOC into disjoint MOCs.
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage in the store.
  /// * ``include_indirect_neighbours`` -
  ///     if `false`, only consider  cells having a common edge as been part of a same MOC
  ///     if `true`, also consider cells having a common vertex as been part of the same MOC
  ///
  /// # Errors
  ///
  /// * ``depth`` is not comprised in `[0, Hpx::<T>::MAX_DEPTH] = [0, 29]`
  #[pyfn(m)]
  fn split(
    index: usize,
    include_indirect_neighbours: bool,
  ) -> PyResult<Vec<usize>> {
    if include_indirect_neighbours {
      U64MocStore::get_global_store().split_indirect(index)
    } else {
      U64MocStore::get_global_store().split(index)
    }.map_err(PyValueError::new_err)
  }

  /// Degrade a (1D) coverage to a specific depth.
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage in the store.
  /// * ``depth`` - The depth to degrade the time coverage to.
  ///
  #[pyfn(m)]
  fn degrade(
    index: usize,
    depth: u8,
  ) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .degrade(index, depth)
      .map_err(PyValueError::new_err)
  }

  #[pyfn(m)]
  fn get_barycenter(py: Python, index: usize) -> PyResult<Py<PyTuple>> {
    U64MocStore::get_global_store()
      .barycenter(index)
      .map(|(lon, lat)| PyTuple::new(py, vec![lon, lat]).into())
      .map_err(PyValueError::new_err)
  }

  #[pyfn(m)]
  fn get_largest_distance_from_coo_to_moc_vertices(index: usize, lon_rad: f64, lat_rad: f64) -> PyResult<f64> {
    U64MocStore::get_global_store()
      .largest_distance_from_coo_to_moc_vertices(index, lon_rad, lat_rad)
      .map_err(PyValueError::new_err)
  }


  /// Get the depth of a spatial coverage.
  ///
  /// # Arguments
  ///
  /// * ``index`` - Index of the coverage in the store.
  #[pyfn(m)]
  fn get_smoc_depth(index: usize) -> PyResult<u8> {
    U64MocStore::get_global_store()
      .get_smoc_depth(index)
      .map_err(PyValueError::new_err)
  }

  /// Get the depth of a time coverage.
  ///
  /// # Arguments
  ///
  /// * ``index`` - Index of the coverage in the store.
  #[pyfn(m)]
  fn get_tmoc_depth(index: usize) -> PyResult<u8> {
    U64MocStore::get_global_store()
      .get_tmoc_depth(index)
      .map_err(PyValueError::new_err)
  }


  /// Compute the sum of all ranges size
  ///
  /// # Arguments
  ///
  /// * ``index`` - Index of the coverage in the store.
  #[pyfn(m)]
  fn ranges_sum(index: usize) -> PyResult<u64> {
    U64MocStore::get_global_store()
      .get_ranges_sum(index)
      .map_err(PyValueError::new_err)
  }

  /// Returns the 1st index in the MOC (error if empty MOC).
  ///
  /// # Arguments
  ///
  /// * ``index`` - Index of the coverage in the store.
  #[pyfn(m)]
  fn first_index(index: usize) -> PyResult<u64> {
    U64MocStore::get_global_store()
      .get_1st_axis_min(index)
      .and_then(|opt| opt.ok_or_else(|| String::from("No min value in an empty MOC")))
      .map_err(PyValueError::new_err)
  }

  /// Returns the last index in the MOC (error if empty MOC).
  ///
  /// # Arguments
  ///
  /// * ``index`` - Index of the coverage in the store.
  #[pyfn(m)]
  fn last_index(index: usize) -> PyResult<u64> {
    U64MocStore::get_global_store()
      .get_1st_axis_max(index)
      .and_then(|opt| opt.ok_or_else(|| String::from("No max value in an empty MOC")))
      .map_err(PyValueError::new_err)
  }


  /// Compute the coverage fraction of a MOC
  ///
  /// # Arguments
  ///
  /// * ``index`` - Index of the coverage in the store.
  #[pyfn(m)]
  fn coverage_fraction(index: usize) -> PyResult<f64> {
    U64MocStore::get_global_store()
      .get_coverage_percentage(index)
      .map(|c| 0.01 * c) // / 100 to transform a percentage in a fraction
      .map_err(PyValueError::new_err)
  }

  /// Get the **uniq** HEALPix cell indices from the MOC of given index.
  ///
  /// # Arguments
  ///
  /// * ``ranges`` - The HEALPix cells defined in the **nested** format.
  #[pyfn(m)]
  fn to_uniq_hpx(py: Python, index: usize) -> PyResult<Py<PyArray1<u64>>> {
    U64MocStore::get_global_store()
      .to_uniq_hpx(index)
      .map(|v| v.into_pyarray(py).to_owned())
      .map_err(PyValueError::new_err)
  }

  /// Get the generic **uniq** cell indices from the MOC of given index.
  ///
  /// # Arguments
  ///
  /// * ``ranges`` - The HEALPix cells defined in the **nested** format.
  #[pyfn(m)]
  fn to_uniq_gen(py: Python, index: usize) -> PyResult<Py<PyArray1<u64>>> {
    U64MocStore::get_global_store()
      .to_uniq_gen(index)
      .map(|v| v.into_pyarray(py).to_owned())
      .map_err(PyValueError::new_err)
  }

  /// Get the zorder **uniq** cell indices from the MOC of given index.
  ///
  /// # Arguments
  ///
  /// * ``ranges`` - The HEALPix cells defined in the **nested** format.
  #[pyfn(m)]
  fn to_uniq_zorder(py: Python, index: usize) -> PyResult<Py<PyArray1<u64>>> {
    U64MocStore::get_global_store()
      .to_uniq_zorder(index)
      .map(|v| v.into_pyarray(py).to_owned())
      .map_err(PyValueError::new_err)
  }

  #[pyfn(m)]
  fn to_ranges(py: Python, index: usize) -> PyResult<Py<PyArray2<u64>>> {
    U64MocStore::get_global_store()
      .to_ranges(index)
      .map_err(PyValueError::new_err)
      .and_then(|mut v| {
        let len = v.len();
        // We could have used Array::from_shape_vec; to be checked: no clone in both cases
        PyArray1::from_vec(py, utils::flatten(&mut v))
          .reshape(Ix2(len, 2_usize))
          .map(|a| a.to_owned())
      })
  }

  #[pyfn(m)]
  fn to_rgba(py: Python, index: usize, size_y: u16) -> PyResult<Py<PyArray3<u8>>> {
    U64MocStore::get_global_store()
      .to_image(index, size_y)
      .map_err(PyValueError::new_err)
      .and_then(|box_u8| {
        PyArray1::from_slice(py, &box_u8)
          .reshape(Ix3(size_y as usize,(size_y << 1) as usize, 4_usize))
          .map(|a| a.to_owned())
        }
      )
  }

  /// Create a temporal coverage from a list of time values expressed in microseconds since
  /// jd origin.
  ///
  /// # Arguments
  ///
  /// * ``depth`` - depth of the MOC
  /// * ``times`` - The list of time values expressed in microseconds since jd=0
  ///
  /// # Errors
  ///
  /// * If the number of ``min_times`` and ``max_times`` do not match.
  #[pyfn(m)]
  fn from_time_in_microsec_since_jd_origin(
    depth: u8,
    times: PyReadonlyArray1<u64>,
  ) -> PyResult<usize> {
    let times = times.as_array().into_iter().cloned();
    U64MocStore::get_global_store()
      .from_microsec_since_jd0(depth, times)
      .map_err(PyValueError::new_err)
  }


  /// Create a temporal coverage from a list of time ranges expressed in jd.
  ///
  /// # Arguments
  ///
  /// * ``depth`` - depth of the MOC
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
  #[pyfn(m)]
  fn from_time_ranges(
    depth: u8,
    min_times: PyReadonlyArray1<f64>,
    max_times: PyReadonlyArray1<f64>,
  ) -> PyResult<usize> {
    let min_times = min_times.as_array().into_iter();
    let max_times = max_times.as_array().into_iter();
    U64MocStore::get_global_store()
      .from_decimal_jd_ranges(depth, min_times.zip(max_times).map(|(min, max)| *min..*max))
      .map_err(PyValueError::new_err)
  }

  /// Create a temporal coverage from a list of time ranges expressed in microseconds since
  /// jd origin.
  ///
  /// # Arguments
  ///
  /// * ``depth`` - depth of the MOC
  /// * ``min_times`` - The list of inf bounds of the time ranges expressed in microseconds since
  ///    jd origin.
  /// * ``max_times`` - The list of sup bounds of the time ranges expressed in microseconds since
  ///    jd origin.
  ///
  /// # Errors
  ///
  /// * If the number of ``min_times`` and ``max_times`` do not match.
  #[pyfn(m)]
  fn from_time_ranges_in_microsec_since_jd_origin(
    depth: u8,
    min_times: PyReadonlyArray1<u64>,
    max_times: PyReadonlyArray1<u64>,
  ) -> PyResult<usize> {
    let min_times = min_times.as_array().into_iter();
    let max_times = max_times.as_array().into_iter();
    U64MocStore::get_global_store()
      .from_microsec_ranges_since_jd0(depth, min_times.zip(max_times).map(|(min, max)| *min..*max))
      .map_err(PyValueError::new_err)
  }

  /// Flatten cells to the moc depth
  ///
  /// # Arguments
  ///
  /// * ``data`` - The spatial coverage
  #[pyfn(m)]
  fn flatten_to_moc_depth(py: Python, index: usize) -> PyResult<Py<PyArray1<u64>>> {
    U64MocStore::get_global_store()
      .flatten_to_moc_depth(index)
      .map(move |v| PyArray1::from_vec(py, v).to_owned())
      .map_err(PyIOError::new_err)
  }
  //

  /// Flatten cells to a specific depth
  ///
  /// # Arguments
  ///
  /// * ``index`` - The index of the coverage.
  #[pyfn(m)]
  fn flatten_to_depth(py: Python, index: usize, depth: u8) -> PyResult<Py<PyArray1<u64>>> {
    U64MocStore::get_global_store()
      .flatten_to_depth(index, depth)
      .map(move |v| PyArray1::from_vec(py, v).to_owned())
    .map_err(PyIOError::new_err)
  }


  #[pyfn(m)]
  fn new_empty_smoc(depth: u8) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .new_empty_smoc(depth)
      .map_err(PyIOError::new_err)
  }

  #[pyfn(m)]
  fn new_empty_tmoc(depth: u8) -> PyResult<usize> {
    U64MocStore::get_global_store()
      .new_empty_tmoc(depth)
      .map_err(PyIOError::new_err)
  }

  /// Create a spatial coverage from a list of HEALPix cell indices.
  ///
  /// # Arguments
  ///
  /// * ``depth`` - The depth of the created MOC
  /// * ``depth`` - The depths of each HEALPix cell indices
  /// * ``pixels`` - A set of HEALPix cell indices
  ///
  /// # Precondition
  ///
  /// ``pixels`` and ``depth`` must be valid. This means that:
  /// * ``depths`` contains values in the range `[0, <T>::MAXDEPTH] = [0, 29]`
  /// * ``pixels`` contains values in the range `[0, 12*4**(depth)]`
  ///
  /// # Errors
  ///
  /// * ``depth`` and ``pixels`` have not the same length.
  #[pyfn(m)]
  fn from_healpix_cells(
    depth: u8,
    depths: PyReadonlyArrayDyn<u8>,
    pixels: PyReadonlyArrayDyn<u64>,
  ) -> PyResult<usize> {
    let depths_shape = depths.shape();
    let pixels_shape = pixels.shape();
    if depths_shape != pixels_shape {
      return Err(PyValueError::new_err(format!("Depths shape different from pixels shape: {:?} != {:?}", depths_shape, pixels_shape)));
    }
    let depths_it = depths.as_array().into_iter().cloned();
    let pixels_it = pixels.as_array().into_iter().cloned();

    U64MocStore::get_global_store()
      .from_hpx_cells(depth, depths_it.zip(pixels_it), None)
      .map_err(PyValueError::new_err)
  }

  Ok(())
}
