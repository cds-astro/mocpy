//! The `coverage` module contains general methods
//! concerning 1d coverages.
//!
//! These are not specific to a particular dimension
//!

use std::fs;
use std::fs::File;
use std::io::{BufWriter, BufReader};
use std::path::Path;
use std::error::Error;
use std::ops::Range;
use std::collections::HashMap;

use ndarray::{Array, Array1, Array2, Axis};

use pyo3::exceptions;
use pyo3::prelude::{Python, PyResult};
use pyo3::types::{PyDict, PyList, PyString};
use pyo3::{PyObject, ToPyObject};

use moc::qty::{MocQty, Hpx, Time};
use moc::elemset::range::{
    MocRanges, HpxRanges, TimeRanges,
    uniq::HpxUniqRanges
};
use moc::moc::{
    RangeMOCIterator, RangeMOCIntoIterator,
    CellMOCIterator, CellMOCIntoIterator,
    CellOrCellRangeMOCIterator, CellOrCellRangeMOCIntoIterator,
    range::{RangeMOC, op::convert::convert_to_u64}
};

use moc::ranges::Ranges;
use moc::deser::fits::{
    ranges_to_fits_ivoa, from_fits_ivoa,
    MocIdxType, MocQtyType, MocType,
    error::FitsError
};
use moc::deser::ascii::from_ascii_ivoa;
use moc::deser::json::from_json_aladin;

use crate::ndarray_fromto::{array2_to_ranges, array2_to_mocranges};


/// Degrade a coverage.
///
/// # Arguments
///
/// * ``coverage`` - The coverage to degrade.
/// * ``depth`` - The depth to degrade the spatial coverage to.
///
pub fn degrade_ranges<Q: MocQty<u64>>(coverage: &mut MocRanges<u64, Q>, depth: u8) -> PyResult<()> {
   coverage.degrade(depth);
   Ok(())
}

/// Flatten cells to a specific depth (depth = max_depth - delta_depth).
///
/// # Arguments
///
/// * ``coverage`` - The coverage
/// * ``dim`` - The dimension of the quantity in the coverage (1 for Time, 2 for Coordinates on the unit sphere).
/// * ``delta_depth`` - Difference between the maximum depth and the target depth
pub fn flatten_hpx_pixels(data: Array2<u64>, depth: u8) -> Array1<u64> {
    let factor = 1 << (2 * (Hpx::<u64>::MAX_DEPTH - depth)) as u64;
    //  let factor = (dim * depth) as u32;

    let flattened_intervals = &data / &Array::from_elem(data.shape(), factor);

    if flattened_intervals.is_empty() {
        // Do not iter on axis 0
        Array::zeros((0,))
    } else {
        // Starting capacity = 2 * number of ranges, but not smaller than 128 (in Java the default is 16)
        let mut flattened_pixels = Vec::<u64>::with_capacity(128.max(data.len() << 1));
        for interval in flattened_intervals.axis_iter(Axis(0)) {
            for pix in interval[0]..interval[1] {
                flattened_pixels.push(pix);
            }
        }
        flattened_pixels.into()
    }
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
pub fn from_json(py: Python, input: &PyDict) -> PyResult<HpxRanges<u64>> {
    const TYPE_KEY_MSG_ERR: &str =
        "The key must be a python str that \n \
         encodes an integer on 1 byte. Ex: {'5': [0, 6, 7, ..., 9]}";
    const TYPE_VALUES_MSG_ERR: &str =
        "The values must be a list of unsigned \n \
         integer on 64 bits (8 bytes). Ex: {'5': [0, 6, 7, ..., 9]}";
    const EXTRACT_IPIX_FROM_LIST_MSG_ERR: &str =
        "Cannot extract 64 bits unsigned integer from the python list!";

    let mut ranges = Vec::<Range<u64>>::new();

    for (depth, pixels) in input {
        let depth = depth
            .downcast::<PyString>()
            .map_err(|_| exceptions::PyTypeError::new_err(TYPE_KEY_MSG_ERR))?
            .to_string()
            .parse::<u8>()
            .unwrap();

        let pixels = pixels
            .downcast::<PyList>()
            .map_err(|_| exceptions::PyTypeError::new_err(TYPE_VALUES_MSG_ERR))?;

        let shift = ((Hpx::<u64>::MAX_DEPTH - depth) << 1) as u64;

        for p in pixels.into_iter() {
            let pixel = p
                .to_object(py)
                .extract::<u64>(py)
                .map_err(|_| exceptions::PyValueError::new_err(EXTRACT_IPIX_FROM_LIST_MSG_ERR))?;

            let e1 = pixel << shift;
            let e2 = (pixel + 1) << shift;
            ranges.push(e1..e2);
        }
    }

    let result = HpxRanges::<u64>::new_from(ranges);
    Ok(result)
}


/// Serializes a spatial coverage to a JSON format
///
/// # Arguments
///
/// * ``coverage`` - The spatial coverage ranges to serialize.
pub fn to_json(py: Python, coverage: HpxRanges<u64>) -> PyResult<PyObject> {
    let result = PyDict::new(py);

    let size = (Hpx::<u64>::MAX_DEPTH + 1) as usize;
    let mut dict = HashMap::with_capacity(size);

    // Initialize the hash map with all the depths and
    // an empty vector of pixels associated to each of them
    for d in 0..size {
        dict.insert(d.to_string(), Vec::<u64>::new()).unwrap();
    }

    // Fill the hash map with the pixels
    for (depth, pix) in coverage.iter_depth_pix() {
        dict.get_mut(&depth.to_string()).unwrap().push(pix);
    }

    // Fill the dictionary with only the depths
    // where pixels are found
    for (d, ipix) in &dict {
        if !ipix.is_empty() {
            result.set_item(d, ipix).map_err(|_| {
                exceptions::PyValueError::new_err(
                    "An error occured when inserting items into the PyDict",
                )
            })?;
        }
    }

    Ok(result.to_object(py))
}

/// Make a coverage consistent.
/// Now, the consistency is ensure before calling this method.
/// We should rename it `divide`, and put the test on `min_depth` in the calling method.
///
/// # Infos
///
/// This is an internal method whose purpose is not to be called
/// by an user. It is called inside of the `mocpy.IntervalSet` class.
///
/// # Arguments
///
/// * ``coverage`` - The coverage to make consistent.
/// * ``min_depth`` - A minimum depth. This argument is optional.
///   A min depth means that we do not want any HEALPix cells to be
///   of depth < to ``min_depth``. This argument is used for example for
///   plotting a spatial coverage. All HEALPix cells of depth < 2 are splitted
///   into cells of depth 2.
///
/// # Errors
///
/// * ``min_depth`` is not comprised in `[0, MocQty::<T>::MAXDEPTH]`
pub fn merge<Q: MocQty<u64>>(mut coverage: MocRanges<u64, Q>, min_depth: i8) -> PyResult<MocRanges<u64, Q>> {
    // If a min depth has been given
    if min_depth != -1 {
        // Then we check its validity
        let max_depth = Q::MAX_DEPTH;
        if min_depth < 0 || min_depth as u8 > max_depth {
            return Err(exceptions::PyValueError::new_err("Min depth is not valid."));
        }

        // And perform the division of the ranges
        coverage = coverage.divide(min_depth as u8);
    }

    Ok(coverage)
}

/// Compute the depth of a coverage
///
/// # Arguments
///
/// * ``coverage`` - The input coverage.
pub fn depth<Q: MocQty<u64>>(coverage: &MocRanges<u64, Q>) -> u8 {
    coverage.compute_min_depth()
}

/// Compute the sky fraction of a spatial coverage
///
/// # Arguments
///
/// * ``coverage`` - The spatial coverage
pub fn sky_fraction(coverage: &Array2<u64>) -> f32 {
    if coverage.is_empty() {
        return 0_f32;
    }

    let sum = coverage.rows()
        .into_iter()
        .fold(0_u64, |acc, r| acc + r[1] - r[0]);

    sum as f32 / ((12_u64 << 58) as f32)
}


/// Cast an `Array2<u64>` coming from MOCPy python code to
/// a generic `Ranges<u64>` object.
pub fn create_ranges_from_py_unchecked(data: Array2<u64>) -> Ranges<u64> {
   array2_to_ranges(data)
}

pub fn build_ranges_from_py(data: Array2<u64>) -> Ranges<u64> {
    Ranges::new_from(create_ranges_from_py_unchecked(data).0.to_vec())
}

/// Cast an `Array2<u64>` coming from MOCPy python code to
/// a `HpxRanges<u64>` object.
pub fn create_hpx_ranges_from_py_unchecked(data: Array2<u64>) -> HpxRanges<u64> {
    array2_to_mocranges(data)
}

pub fn build_hpx_ranges_from_py(data: Array2<u64>) -> HpxRanges<u64> {
    HpxRanges::new_from(create_ranges_from_py_unchecked(data).0.to_vec())
}

/// Cast an `Array2<u64>` coming from MOCPy python code to
/// a `TimeRanges<u64>` object.
pub fn create_time_ranges_from_py_uncheked(data: Array2<u64>) -> TimeRanges<u64> {
    array2_to_mocranges(data)
}

pub fn build_time_ranges_from_py(data: Array2<u64>) -> TimeRanges<u64> {
    TimeRanges::new_from(create_ranges_from_py_unchecked(data).0.to_vec())
}


/// Cast an `Array2<u64>` coming from MOCPy python code to
/// an `HpxUniqRanges<u64>` object.
pub fn create_uniq_ranges_from_py(data: Array2<u64>) -> HpxUniqRanges<u64> {
    HpxUniqRanges::new_from_sorted(create_ranges_from_py_unchecked(data).0.to_vec())
}



// Save/Serialize

pub fn to_ascii_str<Q: MocQty<u64>>(depth_max: u8, ranges: MocRanges<u64, Q>) -> String {
    let moc = RangeMOC::new(depth_max, ranges);
    let mut ascii = Vec::new();
    moc.into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(None,false, &mut ascii)
      .unwrap(); // unwrap since we do not write in a file but in a string
    unsafe{ String::from_utf8_unchecked(ascii) }
}

pub fn to_ascii_file<Q: MocQty<u64>>(depth_max: u8, ranges: MocRanges<u64, Q>, path: String) -> std::io::Result<()> {
    let moc = RangeMOC::new(depth_max, ranges);
    let file = File::create(Path::new(&path))?;
    let writer = BufWriter::new(file);
    moc.into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(Some(80),false, writer)
}

pub fn to_json_str<Q: MocQty<u64>>(depth_max: u8, ranges: MocRanges<u64, Q>) -> String {
    let moc = RangeMOC::new(depth_max, ranges);
    let mut json = Vec::new();
    moc.into_range_moc_iter()
      .cells()
      .to_json_aladin(Some(80), &mut json)
      .unwrap(); // unwrap since we do not write in a file but in a string
    unsafe{ String::from_utf8_unchecked(json) }
}

pub fn to_json_file<Q: MocQty<u64>>(depth_max: u8, ranges: MocRanges<u64, Q>, path: String) -> std::io::Result<()> {
    let moc = RangeMOC::new(depth_max, ranges);
    let file = File::create(Path::new(&path))?;
    let writer = BufWriter::new(file);
    moc.into_range_moc_iter()
      .cells()
      .to_json_aladin(Some(80), writer)
}

pub fn to_fits_file<Q: MocQty<u64>>(depth_max: u8, ranges: MocRanges<u64, Q>, path: String) -> Result<(), FitsError> {
    let moc = RangeMOC::new(depth_max, ranges);
    let file = File::create(Path::new(&path))?;
    let writer = BufWriter::new(file);
    ranges_to_fits_ivoa(moc.into_range_moc_iter(), None, None, writer)
}

// Load/deserialize
pub fn from_fits_file_spatial(path: String) -> Result<MocRanges<u64, Hpx::<u64>>, Box<dyn Error>> {
  let file = File::open(&path).map_err(Box::new)?;
  let reader = BufReader::new(file);
  let ranges: Vec<Range<u64>> = match from_fits_ivoa(reader).map_err(Box::new)? {
    MocIdxType::U16(MocQtyType::Hpx(MocType::Cells(moc))) => convert_to_u64::<_, Hpx<_>, _, Hpx<u64>>(moc.into_cell_moc_iter().ranges()).collect(),
    MocIdxType::U16(MocQtyType::Hpx(MocType::Ranges(moc))) => convert_to_u64::<_, Hpx<_>, _, Hpx<u64>>(moc).collect(),
    MocIdxType::U32(MocQtyType::Hpx(MocType::Cells(moc))) => convert_to_u64::<_, Hpx<_>, _, Hpx<u64>>(moc.into_cell_moc_iter().ranges()).collect(),
    MocIdxType::U32(MocQtyType::Hpx(MocType::Ranges(moc))) => convert_to_u64::<_, Hpx<_>, _, Hpx<u64>>(moc).collect(),
    MocIdxType::U64(MocQtyType::Hpx(MocType::Cells(moc))) => moc.into_cell_moc_iter().ranges().collect(),
    MocIdxType::U64(MocQtyType::Hpx(MocType::Ranges(moc))) => moc.collect(),
    _ => return Err(String::from("FITS file  content not compatible with a space moc of u64").into())
  };
  Ok(MocRanges::new_unchecked(ranges))
}

pub fn from_fits_file_time(path: String) -> Result<MocRanges<u64, Time::<u64>>, Box<dyn Error>> {
    let file = File::open(&path).map_err(Box::new)?;
    let reader = BufReader::new(file);
    let ranges: Vec<Range<u64>> = match from_fits_ivoa(reader).map_err(Box::new)? {
        MocIdxType::U16(MocQtyType::Time(MocType::Cells(moc))) => convert_to_u64::<_, Time<_>, _, Time<u64>>(moc.into_cell_moc_iter().ranges()).collect(),
        MocIdxType::U16(MocQtyType::Time(MocType::Ranges(moc))) => convert_to_u64::<_, Time<_>, _, Time<u64>>(moc).collect(),
        MocIdxType::U32(MocQtyType::Time(MocType::Cells(moc))) => convert_to_u64::<_, Time<_>, _, Time<u64>>(moc.into_cell_moc_iter().ranges()).collect(),
        MocIdxType::U32(MocQtyType::Time(MocType::Ranges(moc))) => convert_to_u64::<_, Time<_>, _, Time<u64>>(moc).collect(),
        MocIdxType::U64(MocQtyType::Time(MocType::Cells(moc))) => moc.into_cell_moc_iter().ranges().collect(),
        MocIdxType::U64(MocQtyType::Time(MocType::Ranges(moc))) => moc.collect(),
        _ => return Err(String::from("FITS file content not compatible with a time moc of u64").into())
    };
    Ok(MocRanges::new_unchecked(ranges))
}

pub fn from_ascii_file<Q: MocQty<u64>>(path: String) -> Result<MocRanges<u64, Q>, Box<dyn Error>> {
    let file_content = fs::read_to_string(&path).map_err(Box::new)?;
    from_ascii_str(file_content)
}

pub fn from_ascii_str<Q: MocQty<u64>>(ascii: String) -> Result<MocRanges<u64, Q>, Box<dyn Error>> {
    let ranges: Vec<Range<u64>> = from_ascii_ivoa::<u64, Q>(&ascii)
      .map_err(Box::new)?
      .into_cellcellrange_moc_iter()
      .ranges()
      .collect();
    Ok(MocRanges::new_unchecked(ranges))
}

pub fn from_json_file<Q: MocQty<u64>>(path: String) -> Result<MocRanges<u64, Q>, Box<dyn Error>> {
    let file_content = fs::read_to_string(&path).map_err(Box::new)?;
    from_json_str(file_content)
}

pub fn from_json_str<Q: MocQty<u64>>(json: String) -> Result<MocRanges<u64, Q>, Box<dyn Error>> {
    let ranges: Vec<Range<u64>> = from_json_aladin::<u64, Q>(&json)?
      .into_cell_moc_iter()
      .ranges()
      .collect();
    Ok(MocRanges::new_unchecked(ranges))
}