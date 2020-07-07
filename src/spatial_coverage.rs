//! The `spatial_coverage` module contains methods
//! related to the creation and manipulation of
//! spatial coverages.

use intervals::nestedranges::NestedRanges;

use intervals::bounded::Bounded;
use std::ops::Range;

use pyo3::exceptions;
use pyo3::prelude::PyResult;
use rayon::prelude::*;

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
    depth: i8,
) -> PyResult<NestedRanges<u64>> {
    if lon.len() != lat.len() {
        return Err(exceptions::ValueError::py_err(
            "Longitudes and Latitudes \
             do not have the same shapes.",
        ));
    }
    if depth < 0 || depth > <u64>::MAXDEPTH {
        return Err(exceptions::ValueError::py_err(format!(
            "Depth must be comprised between in [0, {0}]",
            <u64>::MAXDEPTH
        )));
    }

    let mut data = Vec::<Range<u64>>::with_capacity(lon.len());
    data.resize(lon.len(), 0..1);

    let shift = (<u64>::MAXDEPTH - depth) << 1;
    let layer = healpix::nested::get_or_create(depth as u8);

    data.par_iter_mut()
        .zip_eq(lon.into_par_iter().zip_eq(lat.into_par_iter()))
        .for_each(|(p, (l, b))| {
            let pix = layer.hash(l, b);

            let e1 = pix << shift;
            let e2 = (pix + 1) << shift;
            *p = e1..e2;
        });

    let result = NestedRanges::<u64>::new(data).make_consistent();
    Ok(result)
}

use intervals::valuedcell::valued_cells_to_moc;
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
    max_depth: u32,
    uniq: Array1<u64>,
    values: Array1<f64>,
    cumul_from: f64,
    cumul_to: f64
) -> PyResult<NestedRanges<u64>> {
    if uniq.len() != values.len() {
        Err(
            exceptions::ValueError::py_err(
                "`uniq` and values do not have the same size."
            )
        )
    } else if cumul_from >= cumul_to {
        Err(
            exceptions::ValueError::py_err(
                "`cumul_from` has to be < to `cumul_to`."
            )
        )
    } else {
        // Uniq and values have the same size (can be empty)
        let result = valued_cells_to_moc(max_depth, uniq, values, cumul_from, cumul_to);
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
    coverage_left: &NestedRanges<u64>,
    coverage_right: &NestedRanges<u64>,
) -> NestedRanges<u64> {
    coverage_left.union(coverage_right)
}

/// Performs the intersection between two spatial coverages
///
/// # Arguments
///
/// * `coverage_left` - Left operand
/// * `coverage_right` - Right operand
pub fn intersection(
    coverage_left: &NestedRanges<u64>,
    coverage_right: &NestedRanges<u64>,
) -> NestedRanges<u64> {
    coverage_left.intersection(coverage_right)
}

/// Performs the difference between two spatial coverages
///
/// # Arguments
///
/// * `coverage_left` - Left operand
/// * `coverage_right` - Right operand
pub fn difference(
    coverage_left: &NestedRanges<u64>,
    coverage_right: &NestedRanges<u64>,
) -> NestedRanges<u64> {
    coverage_left.difference(coverage_right)
}

/// Computes the complement od a spatial coverage
///
/// # Arguments
///
/// * `coverage` - The input spatial coverage
pub fn complement(coverage: &NestedRanges<u64>) -> NestedRanges<u64> {
    coverage.complement()
}

use crate::coverage;
use ndarray::{Array1, Array2, Axis, Zip};
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
pub fn from_healpix_cells(mut pixels: Array1<u64>, depth: Array1<i8>) -> PyResult<Array2<u64>> {
    let ones: Array1<u64> = Array1::<u64>::ones(pixels.shape()[0]);
    let mut pixels_1 = &pixels + &ones;

    if pixels.shape() != depth.shape() {
        return Err(exceptions::IndexError::py_err(
            "pixels and depth arrays must have the same shape",
        ));
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

    let ranges = coverage::create_nested_ranges_from_py(ranges).make_consistent();

    let result: Array2<u64> = ranges.into();
    Ok(result)
}

use intervals::uniqranges::UniqRanges;
/// Convert a spatial coverage from the **uniq** to the **nested** format.
///
/// # Arguments
///
/// * ``coverage`` - The spatial coverage defined in the **uniq** format.
pub fn to_nested(coverage: UniqRanges<u64>) -> NestedRanges<u64> {
    coverage.to_nested()
}

/// Convert a spatial coverage from the **nested** to the **uniq** format.
///
/// # Arguments
///
/// * ``coverage`` - The spatial coverage defined in the **nested** format.
pub fn to_uniq(coverage: NestedRanges<u64>) -> UniqRanges<u64> {
    coverage.to_uniq()
}
