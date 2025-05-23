//! Methods releted to SF-MOC

use numpy::{PyArray1, PyArrayMethods, PyReadonlyArrayDyn};
use pyo3::{
  exceptions::{PyIOError, PyValueError},
  prelude::*,
  pyfunction, PyResult,
};

use moc::{qty::Frequency, storage::u64idx::U64MocStore};

/// Computes the depth of a Space-Frequency coverage
///
/// # Arguments
///
/// * ``index`` - The index of the Space-Frequency coverage.
///
/// # Infos
///
/// If the Space-Frequency coverage is empty, the returned
/// depth is `(0, 0)`.
#[pyfunction]
pub fn coverage_sf_depth(index: usize) -> PyResult<(u8, u8)> {
  U64MocStore::get_global_store()
    .get_sfmoc_depths(index)
    .map_err(PyIOError::new_err)
}

#[pyfunction]
pub fn new_empty_sfmoc(depth_frequency: u8, depth_space: u8) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .new_empty_sfmoc(depth_frequency, depth_space)
    .map_err(PyIOError::new_err)
}

/// Returns the minimum frequency value of the Space-Frequency coverage
///
/// # Arguments
///
/// * ``index`` - The index of the Space-Frequency coverage.
///
#[pyfunction]
pub fn coverage_sf_min_freq(index: usize) -> PyResult<f64> {
  U64MocStore::get_global_store()
    .get_1st_axis_min(index)
    .and_then(|opt| opt.ok_or_else(|| String::from("Empty SF-MOC")))
    .map(Frequency::<u64>::hash2freq)
    .map_err(PyValueError::new_err)
}

/// Returns the maximum frequency value of the Space-Frequency coverage
///
/// # Arguments
///
/// * ``index`` - The index of the Space-Frequency coverage.
///
#[pyfunction]
pub fn coverage_sf_max_freq(index: usize) -> PyResult<f64> {
  U64MocStore::get_global_store()
    .get_1st_axis_max(index)
    .and_then(|opt| opt.ok_or_else(|| String::from("Empty SF-MOC")))
    .map(Frequency::<u64>::hash2freq)
    .map_err(PyValueError::new_err)
}

/// Deserialize a Frequency-Space coverage from a FITS file.
///
/// # Arguments
///
/// * ``index`` - the index used to store the Frequency-Space coverage
/// * ``path`` - the FITS file path
///
/// # Warning
///
/// This function is not compatible with pre-v2.0 MOC standard.
///
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in writing the FITS file.
#[pyfunction]
pub fn coverage_sf_from_fits_file(path: String) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .load_sfmoc_from_fits_file(path)
    .map_err(PyIOError::new_err)
}

#[pyfunction]
pub fn sfmoc_from_fits_raw_bytes(raw_bytes: &[u8]) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .load_sfmoc_from_fits_buff(raw_bytes)
    .map_err(PyIOError::new_err)
}

/// Deserialize a Frequency-Space coverage from an ASCII file (compatible with the MOC v2.0 standard).
///
/// # Arguments
///
/// * ``index`` - the index used to store the Frequency-Space coverage
/// * ``path`` - the ASCII file path
///
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in writing the FITS file.
#[pyfunction]
pub fn coverage_sf_from_ascii_file(path: String) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .load_sfmoc_from_ascii_file(path)
    .map_err(PyIOError::new_err)
}

/// Deserialize a Frequency-Space coverage from an JSON file.
///
/// # Arguments
///
/// * ``index`` - the index used to store the Frequency-Space coverage
/// * ``path`` - the JSON file path
///
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in writing the FITS file.
#[pyfunction]
pub fn coverage_sf_from_json_file(path: String) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .load_sfmoc_from_json_file(path)
    .map_err(PyIOError::new_err)
}

/// Deserialize a Frequency-Space coverage from an ASCII string (compatible with the MOC v2.0 standard).
///
/// # Arguments
///
/// * ``index`` - the index used to store the Frequency-Space coverage
/// * ``ascii`` - the ASCII string
///
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in writing the FITS file.
#[pyfunction]
pub fn coverage_sf_from_ascii_str(ascii: String) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .load_sfmoc_from_ascii(&ascii)
    .map_err(PyIOError::new_err)
}

/// Deserialize a Frequency-Space coverage from an JSON string.
///
/// # Arguments
///
/// * ``index`` - the index used to store the Frequency-Space coverage
/// * ``json`` - the JSON string
///
/// # Errors
///
/// This method returns a `PyIOError` if the the function fails in writing the FITS file.
#[pyfunction]
pub fn coverage_sf_from_json_str(json: String) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .load_sfmoc_from_json(&json)
    .map_err(PyIOError::new_err)
}

/// Create a 2D Frequency-Space coverage from a list of
/// (frequency, longitude, latitude) tuples.
///
/// # Arguments
///
/// * ``frequencies_hz`` - The frequencies associated to the sky coordinates, in Hz.
/// * ``freq_depth`` - The depth along the Frequency axis.
/// * ``lon_rad`` - The longitudes in radians
/// * ``lat_rad`` - The latitudes in radians
/// * ``pos_depth`` - The depth along the Space axis.
///
/// # Precondition
///
/// * ``lon_rad`` and ``lat_rad`` must be expressed in radians.
/// * ``frequencies_hz`` must be expressed in Hz.
///
/// # Errors
///
/// * ``lon_rad``, ``lat_rad`` and ``Frequencies_hz`` do not have the same length.
/// * ``freq_depth`` is not comprised in `[0, <F>::MAXDEPTH] = [0, 59]`
/// * ``space_depth`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
#[pyfunction]
pub fn from_freq_lonlat(
  frequencies_hz: PyReadonlyArrayDyn<f64>,
  freq_depth: u8,
  lon_rad: PyReadonlyArrayDyn<f64>,
  lat_rad: PyReadonlyArrayDyn<f64>,
  pos_depth: u8,
) -> PyResult<usize> {
  let frequencies_hz = frequencies_hz.to_vec().map_err(PyValueError::new_err)?;
  let lon_rad = lon_rad.to_vec().map_err(PyValueError::new_err)?;
  let lat_rad = lat_rad.to_vec().map_err(PyValueError::new_err)?;

  U64MocStore::get_global_store()
    .create_from_hz_positions(frequencies_hz, lon_rad, lat_rad, freq_depth, pos_depth)
    .map_err(PyValueError::new_err)
}

/// Create a 2D Frequency-Space coverage from a list of
/// (frequency_range, longitude, latitude) tuples.
///
/// # Arguments
///
/// * ``freq_hz_min`` - The beginning Frequency of observation, in Hz.
/// * ``freq_hz_max`` - The ending Frequency of observation, in Hz.
/// * ``freq_depth`` - The depth along the Frequency axis.
/// * ``lon_rad`` - The longitudes in radians
/// * ``lat_rad`` - The latitudes in radians
/// * ``hpx_depth`` - The depth along the Space axis.
///
/// # Precondition
///
/// * ``lon_rad`` and ``lat_rad`` must be expressed in radians.
/// * ``freq_hz_min`` and ``freq_hz_max`` must be expressed in Hz.
///
/// # Errors
///
/// * ``lon_rad``, ``lat_rad``, ``freq_hz_min`` and ``freq_hz_max`` do not have the same length.
/// * ``freq_depth`` is not comprised in `[0, <F>::MAXDEPTH] = [0, 59]`
/// * ``hpx_depth`` is not comprised in `[0, <S>::MAXDEPTH] = [0, 29]`
///
#[pyfunction]
pub fn from_freq_ranges_lonlat(
  freq_hz_min: PyReadonlyArrayDyn<f64>,
  freq_hz_max: PyReadonlyArrayDyn<f64>,
  freq_depth: u8,
  lon_rad: PyReadonlyArrayDyn<f64>,
  lat_rad: PyReadonlyArrayDyn<f64>,
  hpx_depth: u8,
) -> PyResult<usize> {
  let freq_hz_min = freq_hz_min.to_vec().map_err(PyValueError::new_err)?;
  let freq_hz_max = freq_hz_max.to_vec().map_err(PyValueError::new_err)?;
  let lon_rad = lon_rad.to_vec().map_err(PyValueError::new_err)?;
  let lat_rad = lat_rad.to_vec().map_err(PyValueError::new_err)?;

  U64MocStore::get_global_store()
    .create_from_hzranges_positions(
      freq_hz_min,
      freq_hz_max,
      lon_rad,
      lat_rad,
      freq_depth,
      hpx_depth,
    )
    .map_err(PyValueError::new_err)
}

/// Create a Space-Frequency coverage from a list of
/// (frequency_range, longitude, latitude, radius) tuples.
///
/// # Arguments
///
/// * ``frequencies_min`` - Minimum value for the frequency, in Hz
/// * ``frequencies_max`` - Maximum value for the frequency, in Hz
/// * ``d1`` - The depth along the frequency axis.
/// * ``spatial_coverages`` - List of spatial coverages.
///
#[pyfunction]
pub fn from_frequency_ranges_spatial_coverages(
  frequencies_min: PyReadonlyArrayDyn<f64>,
  frequencies_max: PyReadonlyArrayDyn<f64>,
  d1: u8,
  spatial_coverages: PyReadonlyArrayDyn<usize>,
) -> PyResult<usize> {
  let frequencies_min = frequencies_min.to_vec().map_err(PyValueError::new_err)?;
  let frequencies_max = frequencies_max.to_vec().map_err(PyValueError::new_err)?;
  let spatial_coverage_indices = spatial_coverages.to_vec().map_err(PyValueError::new_err)?;
  U64MocStore::get_global_store()
    .from_freq_ranges_spatial_coverages_in_store(
      frequencies_min,
      frequencies_max,
      d1,
      spatial_coverage_indices,
    )
    .map_err(PyValueError::new_err)
}

/// Check if (time, position) tuples are contained into a Time-Space coverage
///
/// # Arguments
///
/// * ``index`` - The index of the Time-Space coverage.
/// * ``frequencies`` - Frequencies in Hz
/// * ``lon`` - The longitudes in radians.
/// * ``lat`` - The latitudes in radians.
///
/// # Errors
///
/// * If `lon`, `lat` and `times` do not have the same length
#[pyfunction]
#[pyo3(pass_module)]
pub fn sfmoc_contains<'a, 'py>(
  module: &Bound<'py, PyModule>,
  index: usize,
  frequencies: PyReadonlyArrayDyn<'a, f64>,
  lon: PyReadonlyArrayDyn<'a, f64>,
  lat: PyReadonlyArrayDyn<'a, f64>,
) -> PyResult<Bound<'py, PyArray1<bool>>> {
  let it_freq = frequencies.as_array().into_iter().cloned();
  let it_lon = lon.as_array().into_iter().cloned();
  let it_lat = lat.as_array().into_iter().cloned();
  let it = it_freq.zip(it_lon.zip(it_lat));
  U64MocStore::get_global_store()
    .filter_freqpos(index, it, |b| b)
    .map(|vec_bool| PyArray1::<bool>::from_vec(module.py(), vec_bool))
    .map_err(PyIOError::new_err)
}

/// Project the Frequency-Space coverage into its second dimension
/// (i.e. the Space axis)
///
/// # Arguments
///
/// * ``fmoc_index`` - The index of the Frequency coverage.
/// * ``sfmoc_index`` - The index of the Frequency-Space coverage.
///
/// # Algorithm
///
/// Returns the union of the spatial coverages for which
/// their frequency ranges are contained into ``x``.
///
/// # Panic
///
/// If the ``ranges`` is not valid i.e.:
///
/// * Contains ranges for which their inf bound is
///   superior to their sup bound.
///
/// This **should** not panic as this code is wrapped around MOCPy
#[pyfunction]
pub fn project_on_sfmoc_space_dim(fmoc_index: usize, sfmoc_index: usize) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .frequency_fold(fmoc_index, sfmoc_index)
    .map_err(PyValueError::new_err)
}

/// Returns the index of a F-MOC
#[pyfunction]
pub fn project_on_sfmoc_freq_dim(smoc_index: usize, sfmoc_index: usize) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .space_fold(smoc_index, sfmoc_index)
    .map_err(PyValueError::new_err)
}
