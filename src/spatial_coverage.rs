//! The `spatial_coverage` module contains methods
//! related to the creation and manipulation of
//! spatial coverages.

use std::ops::Range;

use rayon::prelude::*;

use pyo3::exceptions;
use pyo3::prelude::PyResult;

use moc::qty::{MocQty, Hpx};
use moc::elem::valuedcell::{valued_cells_to_moc, valued_cells_to_moc_with_opt};
use moc::elemset::range::{
    HpxRanges,
    uniq::HpxUniqRanges
};
use moc::ranges::SNORanges;
use moc::moc::range::RangeMOC;
use moc::deser::fits::error::FitsError;

use super::coverage;
use super::ndarray_fromto::mocranges_to_array2;
use ndarray::{ArrayD, ArrayViewD};

use ndarray::{Array1, Array2, Zip};

/*
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
pub fn from_healpix_map(depth: u8, mut pixels: Array1<u64>) -> PyResult<Array2<u64>> {
    let factor = (Hpx::<u64>::MAX_DEPTH - depth ) << 1;
    let ones: Array1<u64> = Array1::<u64>::ones(pixels.shape()[0]);
    let mut pixels_1 = &pixels + &ones; // Probably better to clone and +1 in par_map_inplace
    let to_max_depth = |pix: &mut u64| *pix <<= factor;
    pixels.par_map_inplace(to_max_depth);
    pixels_1.par_map_inplace(to_max_depth);
    Ok(from_lower_and_upperd_bounds(pixels, pixels_1))
}

fn from_lower_and_upperd_bounds(low: Array1<u64>, upp: Array1<u64>) -> Array2<u64> {
    let shape = (low.shape()[0], 1);
    let low = low.into_shape(shape).unwrap();
    let upp = upp.into_shape(shape).unwrap();
    debug_assert_eq!(low.len(), upp.len());
    let mut ranges: Vec<Range<u64>> = Vec::with_capacity(low.len());
    for (start, end) in low.into_iter().zip(upp.into_iter()) {
        ranges.push(start..end);
    }
    mocranges_to_array2(HpxRanges::<u64>::new_from(ranges))
}
*/

/// Creates a spatial coverage from a list of sky coordinates.
///
/// # Arguments
///
/// * ``lon`` - The longitudes of the sky coordinates.
/// * ``lat`` - The latitudes of the sky coordinates.
/// * ``ds`` - The depth at which HEALPix cell indices
///   will be computed.
///
/// # Precondition
///
/// * ``lon`` and ``lat`` are expressed in radians.
/// They are valid because they come from
/// `astropy.units.Quantity` objects.
///
/// # Errors
///
/// If the number of longitudes and latitudes do not match.
pub fn contains(
    coverage: &HpxRanges<u64>,
    lon: ArrayViewD<f64>,
    lat: ArrayViewD<f64>,
    result: &mut ArrayD<bool>,
) -> PyResult<()> {
    // Retrieve the spatial depth of the Space coverage
    let layer = healpix::nested::get(Hpx::<u64>::MAX_DEPTH);
    Zip::from(result)
      .and(&lon)
      .and(&lat)
      .par_for_each(|r, &l, &b| {
          // Compute the HEALPix cell range at the max depth
          // along the spatial dimension
          let pix = layer.hash(l, b);
          
          // Check whether the (time in µs, HEALPix cell nested range)
          // is contained into the Time coverage
          *r = coverage.contains_val(&pix);
      });

    Ok(())
}


/// Convert a spatial coverage from the **uniq** to the **nested** format.
///
/// # Arguments
///
/// * ``coverage`` - The spatial coverage defined in the **uniq** format.
pub fn to_nested(coverage: HpxUniqRanges<u64>) -> HpxRanges<u64> {
    coverage.into_hpx()
}
