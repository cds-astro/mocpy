//! Methods releted to SF-MOC

use numpy::{PyArrayMethods, PyReadonlyArrayDyn};
use pyo3::{
  exceptions::{PyIOError, PyValueError},
  pyfunction, PyResult,
};

use moc::storage::u64idx::U64MocStore;

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
    .time_fold(fmoc_index, sfmoc_index)
    .map_err(PyValueError::new_err)
}

/// Returns the index of a F-MOC
#[pyfunction]
pub fn project_on_sfmoc_freq_dim(smoc_index: usize, sfmoc_index: usize) -> PyResult<usize> {
  U64MocStore::get_global_store()
    .space_fold(smoc_index, sfmoc_index)
    .map_err(PyValueError::new_err)
}
