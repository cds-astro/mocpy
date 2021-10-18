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

/// Create a spatial coverage from a list of sky coordinates
///
/// # Arguments
///
/// * ``lon`` - The longitudes of the sky coordinates.
/// * ``lat`` - The latitudes of the sky coordinates.
/// * ``depth`` - The depth at which HEALPix cell indices
///   will be computed. This will correspond to the depth
///   of the coverage once it will be created.
///
/// # Precondition
///
/// ``lon`` and ``lat`` are expressed in radians.
/// They are valid because they come from
/// `astropy.units.Quantity` objects.
///
/// # Errors
///
/// * If the number of longitudes and latitudes do not match.
/// * If ``depth`` is not in `[0, <u64>::MAXDEPTH] = [0, 29]`
pub fn create_from_position(
    lon: Vec<f64>,
    lat: Vec<f64>,
    depth: u8,
) -> PyResult<HpxRanges<u64>> {
    if lon.len() != lat.len() {
        return Err(exceptions::PyValueError::new_err(
            "Longitudes and Latitudes \
             do not have the same shapes.",
        ));
    }
    if depth > Hpx::<u64>::MAX_DEPTH {
        return Err(exceptions::PyValueError::new_err(format!(
            "Depth must be comprised between in [0, {0}]",
            Hpx::<u64>::MAX_DEPTH
        )));
    }

    let mut data = Vec::<Range<u64>>::with_capacity(lon.len());
    data.resize(lon.len(), 0..1);

    let shift = (Hpx::<u64>::MAX_DEPTH - depth) << 1;
    let layer = healpix::nested::get(depth as u8);

    data.par_iter_mut()
        .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
        .for_each(|(p, (l, b))| {
            let pix = layer.hash(l, b);

            let e1 = pix << shift;
            let e2 = (pix + 1) << shift;
            *p = e1..e2;
        });

    let result = HpxRanges::<u64>::new_from(data);
    Ok(result)
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
/// * ``max_depth`` - the largest depth of the output MOC, which must be larger or equals to the largest
///   depth in the `uniq` values
///
/// # Precondition
///
/// * ``uniq`` and ``values`` must be of the same size
pub fn from_valued_healpix_cells(
    max_depth: u8,
    uniq: Array1<u64>,
    values: Array1<f64>,
    cumul_from: f64,
    cumul_to: f64
) -> PyResult<HpxRanges<u64>> {
    if uniq.len() != values.len() {
        Err(
            exceptions::PyValueError::new_err(
                "`uniq` and values do not have the same size."
            )
        )
    } else if cumul_from >= cumul_to {
        Err(
            exceptions::PyValueError::new_err(
                "`cumul_from` has to be < to `cumul_to`."
            )
        )
    } else {
        // Uniq and values have the same size (can be empty)
        let result = valued_cells_to_moc(max_depth, uniq.iter(), values.iter(), cumul_from, cumul_to);
        Ok(result)
    }
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
/// * ``max_depth`` - the largest depth of the output MOC, which must be larger or equals to the largest
///   depth in the `uniq` values
/// * `asc`: cumulative value computed from lower to highest densities instead of from highest to lowest
/// * `strict`: (sub-)cells overlapping the `cumul_from` or `cumul_to` values are not added
/// * `no_split`: cells overlapping the `cumul_from` or `cumul_to` values are not recursively split
/// * `reverse_decent`: perform the recursive decent from the highest cell number to the lowest (to be compatible with Aladin)
/// 
/// # Precondition
///
/// * ``uniq`` and ``values`` must be of the same size
pub fn from_valued_healpix_cells_with_opt(
    max_depth: u8,
    uniq: Array1<u64>,
    values: Array1<f64>,
    cumul_from: f64,
    cumul_to: f64,
    asc: bool,
    strict: bool,
    no_split: bool,
    reverse_decent: bool,
) -> PyResult<HpxRanges<u64>> {
    use std::f64::consts::PI;
    if uniq.len() != values.len() {
        Err(
            exceptions::PyValueError::new_err(
                "`uniq` and values do not have the same size."
            )
        )
    } else if cumul_from >= cumul_to {
        Err(
            exceptions::PyValueError::new_err(
                "`cumul_from` has to be < to `cumul_to`."
            )
        )
    } else {
        // Uniq and values have the same size (can be empty)
        let area_per_cell = (PI / 3.0) / (1_u64 << (max_depth << 1) as u32) as f64;  // = 4pi / (12*4^depth)
        let uniq_val_dens: Vec<(u64, f64, f64)> = uniq.iter().zip(values.iter())
          .map(|(uniq, val)| {
              let (cdepth, _) = Hpx::<u64>::from_uniq_hpx(*uniq);
              let n_sub_cells = (1_u64 << (((max_depth - cdepth) << 1) as u32)) as f64;
              (*uniq, *val, val / (n_sub_cells * area_per_cell))
          }).collect();
        
        let result = valued_cells_to_moc_with_opt(
            max_depth,
            uniq_val_dens,
            cumul_from,
            cumul_to,
            asc,
            strict,
            no_split,
            reverse_decent,
        );
        Ok(result)
    }
}

/// Performs the union between two spatial coverages
///
/// # Arguments
///
/// * `coverage_left` - Left operand
/// * `coverage_right` - Right operand
pub fn union(
    coverage_left: &HpxRanges<u64>,
    coverage_right: &HpxRanges<u64>,
) -> HpxRanges<u64> {
    coverage_left.union(coverage_right)
}

/// Performs the intersection between two spatial coverages
///
/// # Arguments
///
/// * `coverage_left` - Left operand
/// * `coverage_right` - Right operand
pub fn intersection(
    coverage_left: &HpxRanges<u64>,
    coverage_right: &HpxRanges<u64>,
) -> HpxRanges<u64> {
    coverage_left.intersection(coverage_right)
}

/// Performs the difference between two spatial coverages
///
/// # Arguments
///
/// * `coverage_left` - Left operand
/// * `coverage_right` - Right operand
pub fn difference(
    coverage_left: &HpxRanges<u64>,
    coverage_right: &HpxRanges<u64>,
) -> HpxRanges<u64> {
    coverage_left.difference(coverage_right)
}

/// Computes the complement od a spatial coverage
///
/// # Arguments
///
/// * `coverage` - The input spatial coverage
pub fn complement(coverage: &HpxRanges<u64>) -> HpxRanges<u64> {
    coverage.complement()
}


/// Extend the given MOC adding an edge of cells at max order.
///
/// # Arguments
/// 
/// * `depth_max` - MOC depth
/// * `coverage` - The input spatial coverage
pub fn expand(depth_max: u8, coverage: HpxRanges<u64>) -> HpxRanges<u64> {
    let moc = RangeMOC::<u64, Hpx<u64>>::new(depth_max, coverage);
    let moc_out = moc.expanded();
    moc_out.into_moc_ranges()
}

/// Contract the given MOC adding an edge of cells at max order.
///
/// # Arguments
/// 
/// * `depth_max` - MOC depth
/// * `coverage` - The input spatial coverage
pub fn contract(depth_max: u8, coverage: HpxRanges<u64>) -> HpxRanges<u64> {
    let moc = RangeMOC::<u64, Hpx<u64>>::new(depth_max, coverage);
    moc.contracted().into_moc_ranges()
}

use ndarray::{Array1, Array2, Zip};
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
pub fn from_healpix_cells(mut pixels: Array1<u64>, depth: Array1<u8>) -> PyResult<Array2<u64>> {
    let ones: Array1<u64> = Array1::<u64>::ones(pixels.shape()[0]);
    let mut pixels_1 = &pixels + &ones;

    if pixels.shape() != depth.shape() {
        return Err(exceptions::PyIndexError::new_err(
            "pixels and depth arrays must have the same shape",
        ));
    }

    // ndarray 15.2: par_apply -> par_for_each
    //   see https://docs.rs/ndarray/0.15.2/ndarray/struct.Zip.html#method.par_for_each
    Zip::from(&mut pixels)
        .and(&mut pixels_1)
        .and(&depth)
        .par_for_each(|pix, pix1, &d| {
            let factor = (Hpx::<u64>::MAX_DEPTH - d) << 1;
            *pix <<= factor;
            *pix1 <<= factor;
        });
    Ok(from_lower_and_upperd_bounds(pixels, pixels_1))
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


/// Convert a spatial coverage from the **uniq** to the **nested** format.
///
/// # Arguments
///
/// * ``coverage`` - The spatial coverage defined in the **uniq** format.
pub fn to_nested(coverage: HpxUniqRanges<u64>) -> HpxRanges<u64> {
    coverage.into_hpx()
}

/// Convert a spatial coverage from the **nested** to the **uniq** format.
///
/// # Arguments
///
/// * ``coverage`` - The spatial coverage defined in the **nested** format.
pub fn to_uniq(coverage: HpxRanges<u64>) -> HpxUniqRanges<u64> {
    coverage.into_hpx_uniq()
}

pub fn to_ascii_str(depth_max: u8, ranges: HpxRanges<u64>) -> String {
    coverage::to_ascii_str(depth_max, ranges)
}
pub fn to_ascii_file(depth_max: u8, ranges: HpxRanges<u64>, path: String) -> std::io::Result<()> {
    coverage::to_ascii_file(depth_max, ranges, path)
}
pub fn to_json_str(depth_max: u8, ranges: HpxRanges<u64>) -> String {
    coverage::to_json_str(depth_max, ranges)
}
pub fn to_json_file(depth_max: u8, ranges: HpxRanges<u64>, path: String) -> std::io::Result<()> {
    coverage::to_json_file(depth_max, ranges, path)
}
pub fn to_fits_file(depth_max: u8, ranges: HpxRanges<u64>, path: String) -> Result<(), FitsError> {
    coverage::to_fits_file(depth_max, ranges, path)
}

pub fn from_ascii_str(ascii: String) -> PyResult<HpxRanges<u64>> {
    coverage::from_ascii_str(ascii).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_ascii_file(path: String) -> PyResult<HpxRanges<u64>> {
    coverage::from_ascii_file(path).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_json_str(json: String) -> PyResult<HpxRanges<u64>> {
    coverage::from_json_str(json).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_json_file(path: String) -> PyResult<HpxRanges<u64>> {
    coverage::from_json_file(path).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
pub fn from_fits_file(path: String) -> PyResult<HpxRanges<u64>> {
    coverage::from_fits_file_spatial(path).map_err(|e| exceptions::PyIOError::new_err(e.to_string()))
}
