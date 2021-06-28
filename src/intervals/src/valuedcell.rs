
use ndarray::Array1;

use std::ops::Range;
use std::cmp::Ordering::Equal;

use crate::ranges::Idx;
use crate::mocqty::Hpx;
use crate::mocrange::HpxRange;


/*
/// Creates a MOC from the given list of uniq cells numbers according to the value they contains.
/// We assume that the value is directly proportional to the covered area (like a flux or a probability).
/// Limits are put to select an area having a cumulative value ranging from a given lower limit
/// to a given upper limit.
/// An example is the selection of a region having between 10 and 90 percent of a flux, or
/// an 90 percent completeness.
///
/// # Precondition
/// * `uniq` and `values` do have the same size.
/// * `uniq` and `values` are not empty.
/// * `cumul_from` < `cumul_to`
///
/// # Errors
/// * if `max_depth` is not > to the finest depth found in the `uniq` cells.
///
/// # Args
/// * `max_depth`: the largest depth of the output MOC, which must be larger or equals to the largest
/// depth in the `uniq` values
/// * `uniq`: the list of uniq cells (i.e. values encoding both the HEALPix depth and cell number)
/// * `values`: values associated to each uniq.
/// * `cumul_from`: the cumulative value from which cells are put in the MOC
/// * `cumul_to`: the cumulative value to which cells are put in the MOC
pub fn valued_cells_to_moc(
    max_depth: u8,
    uniq: Array1<u64>,
    values: Array1<f64>,
    cumul_from: f64,
    cumul_to: f64
) -> HpxRanges<u64> {
    let mut valued_uniq_sorted: Vec<(u64, f64, f64)> = uniq.iter().zip(values.iter())
        .map(|(uniq, val)| {
            let (depth, _icell) = Hpx::<u64>::from_uniq_hpx(*uniq);
            let n_sub_cells = 1_u64 << ((max_depth - depth) << 1);
            (*uniq, *val, val / (n_sub_cells as f64))
        })
    .collect::<Vec<(u64, f64, f64)>>();
    // We use b.comp(a) instead of a.cmp(b) to get the DESC order
    valued_uniq_sorted.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Equal));
    let mut result: Vec<Range<u64>> = Vec::with_capacity(valued_uniq_sorted.len());

    let mut i = 0_usize;
    let mut acc = 0_f64;
    while i < valued_uniq_sorted.len() && acc + valued_uniq_sorted[i].1 <= cumul_from {
        acc = acc + valued_uniq_sorted[i].1;
        result.push(Hpx::<u64>::uniq_hpx_to_range(valued_uniq_sorted[i].0));
        i += 1;
    }
    if i < valued_uniq_sorted.len() && acc < cumul_from {
        let (depth, icell) = Hpx::<u64>::from_uniq_hpx(valued_uniq_sorted[i].0);
        result = recursive_descent(
            depth, icell, max_depth,
            valued_uniq_sorted[i].1, cumul_from - acc, result);
        i += 1;
    }
    while i < valued_uniq_sorted.len() && acc + valued_uniq_sorted[i].1 <= cumul_to {
        acc = acc + valued_uniq_sorted[i].1;
        result.push(Hpx::<u64>::uniq_hpx_to_range(valued_uniq_sorted[i].0));
        i += 1;
    }
    if i < valued_uniq_sorted.len() && acc < cumul_to {
        let (depth, icell) = Hpx::<u64>::from_uniq_hpx(valued_uniq_sorted[i].0);
        result = recursive_descent(
            depth, icell, max_depth,
            valued_uniq_sorted[i].1, cumul_to - acc, result);
    }
    HpxRanges::new(result).make_consistent()
}

fn recursive_descent(
    depth: u8,
    ipix: u64,
    max_depth: u8,

    cell_val: f64,
    mut target_val: f64,

    mut result: Vec<Range<u64>>
) -> Vec<Range<u64>> {
    if depth == max_depth {
        if cell_val <= target_val {
            let rng: HpxRange<u64> = (depth as u8, ipix).into();
            result.push(rng.0);
        }
    } else if target_val > 0_f64 {
        let subcell_val = cell_val / 4_f64;
        let depth = depth + 1;
        let ipix = ipix << 2;
        let mut i = 0_u64;
        while i < 4 && target_val - subcell_val >= 0_f64 {
            let rng: HpxRange<u64> = (depth as u8, ipix + i).into();
            result.push(rng.0);
            target_val = target_val - subcell_val;
            i = i + 1;
        }
        if i < 4 {
            result = recursive_descent(
                depth, ipix + i, max_depth,
                subcell_val, target_val, result
            );
        }
    }
    result
}
*/


// Version using generics. Seem to not work on windows archs. Is it a PyO3 bug or a
// bug involving the windows image used by AppVeyor?

use num::{Num, One};
use crate::mocranges::HpxRanges;

pub trait DivBy<Rhs = Self> {
    type Output;
    fn div_by(self, rhs: Rhs) -> Self::Output;
}

impl DivBy<u64> for f64 {
    type Output = Self;

    fn div_by(self, rhs: u64) -> Self::Output {
        self / (rhs as f64)
    }
}

/// Creates a MOC from the given list of uniq cells numbers according to the value they contains.
/// We assume that the value is directly proportional to the covered area (like a flux or a probability).
/// Limits are put to select an area having a cumulative value ranging from a given lower limit
/// to a given upper limit.
/// An example is the selection of a region having between 10 and 90 percent of a flux, or
/// an 90 percent completeness.
///
/// # Precondition
/// * `uniq` and `values` do have the same size.
/// * `uniq` and `values` are not empty.
/// * `cumul_from` < `cumul_to`
///
/// # Errors
/// * if `max_depth` is not > to the finest depth found in the `uniq` cells.
///
/// # Args
/// * `max_depth`: the largest depth of the output MOC, which must be larger or equals to the largest
/// depth in the `uniq` values
/// * `uniq`: the list of uniq cells (i.e. values encoding both the HEALPix depth and cell number)
/// * `values`: values associated to each uniq.
/// * `cumul_from`: the cumulative value from which cells are put in the MOC
/// * `cumul_to`: the cumulative value to which cells are put in the MOC
pub fn valued_cells_to_moc<T, V>(
    max_depth: u8,
    uniq: Array1<T>,
    values: Array1<V>,
    cumul_from: V,
    cumul_to: V
    ) -> HpxRanges<T>
where T: Idx,
      V: Num + PartialOrd + DivBy<T, Output=V> + Copy + Send + Sync + std::fmt::Debug {
    let mut valued_uniq_sorted: Vec<(T, V, V)> = uniq.iter().zip(values.iter())
        .map(|(uniq, val)| {
            let (depth, _icell) = Hpx::<T>::from_uniq_hpx(*uniq);
            let n_sub_cells = T::one().unsigned_shl(((max_depth - depth) << 1) as u32);
            (*uniq, *val, val.div_by(n_sub_cells))
        })
    .collect::<Vec<(T, V, V)>>();
    // We use b.comp(a) instead of a.cmp(b) to get the DESC order


    valued_uniq_sorted.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Equal));
    let mut result: Vec<Range<T>> = Vec::with_capacity(valued_uniq_sorted.len());

    let mut i = 0_usize;
    let mut acc = V::zero();
    while i < valued_uniq_sorted.len() && acc.add(valued_uniq_sorted[i].1) <= cumul_from {
        acc = acc.add(valued_uniq_sorted[i].1);
        result.push(Hpx::<T>::uniq_hpx_to_range(valued_uniq_sorted[i].0));
        i += 1;
    }

    if i < valued_uniq_sorted.len() && acc < cumul_from {
        let (depth, icell) = Hpx::<T>::from_uniq_hpx(valued_uniq_sorted[i].0);
        result = recursive_descent(
            depth, icell, max_depth,
            valued_uniq_sorted[i].1, cumul_from.sub(acc), result);
        i += 1;
    }
    while i < valued_uniq_sorted.len() && acc.add(valued_uniq_sorted[i].1) <= cumul_to {
        acc = acc.add(valued_uniq_sorted[i].1);
        result.push(Hpx::<T>::uniq_hpx_to_range(valued_uniq_sorted[i].0));
        i += 1;
    }
    if i < valued_uniq_sorted.len() && acc < cumul_to {
        let (depth, icell) = Hpx::<T>::from_uniq_hpx(valued_uniq_sorted[i].0);
        result = recursive_descent(
            depth, icell, max_depth,
            valued_uniq_sorted[i].1, cumul_to.sub(acc), result);
    }
    HpxRanges::new_from_sorted(result)
}

fn recursive_descent<T, V>(
    depth: u8,
    ipix: T,
    max_depth: u8,
    cell_val: V,
    mut target_val: V,
    mut result: Vec<Range<T>>
    ) -> Vec<Range<T>>
where T: Idx,
      V: Num + PartialOrd + DivBy<T, Output=V> + Copy + Send + Sync {
    if depth == max_depth {
        if cell_val <= target_val {
            let rng: HpxRange<T> = (depth, ipix).into();
            result.push(rng.0);
        }
    } else if target_val > V::zero() {
        let four = T::one().unsigned_shl(2);
        let subcell_val = cell_val.div_by(four);
        let depth = depth + 1;
        let ipix = ipix << 2;
        let mut i = T::zero();
        while i < four && target_val.sub(subcell_val) >= V::zero() {
            let rng: HpxRange<T> = (depth, ipix + i).into();
            result.push(rng.0);
            target_val = target_val.sub(subcell_val);
            i = i + One::one();
        }
        if i < four {
            result = recursive_descent(
                depth, ipix + i, max_depth,
                subcell_val, target_val, result
            );
        }
    }
    result
}


#[cfg(test)]
mod tests {
    use std::u64;
    use ndarray::prelude::*;

    use crate::mocqty::{MocQty, Hpx};
    use crate::mocranges::HpxRanges;
    use crate::valuedcell::valued_cells_to_moc;

    #[test]
    fn test_single_uniq() {
        let uniq = array![4];
        let values = array![1_f64];

        let max_depth = 2;

        // let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 0.25_f64);
        let nested_ranges = valued_cells_to_moc(max_depth, uniq, values, 0_f64, 0.25_f64);

        let tdd = ((Hpx::<u64>::MAX_DEPTH - max_depth) << 1) as u32;
        let expect_nested_ranges = HpxRanges::new_unchecked(
            vec![
                0..(4 << tdd)
            ]
        );

        assert_eq!(nested_ranges, expect_nested_ranges);
    }

    #[test]
    fn test_empty() {
        let uniq = array![];
        let values = array![];

        let max_depth = 2;

        // let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 1_f64);
        let nested_ranges = valued_cells_to_moc(max_depth, uniq, values, 0_f64, 1_f64);
        let expect_nested_ranges = HpxRanges::default();

        assert_eq!(nested_ranges, expect_nested_ranges);
    }

    #[test]
    fn test_full_space() {
        let uniq = array![4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        let values = array![0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0_f64, 0_f64];

        let max_depth = 2;

        // let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 1_f64);
        let nested_ranges = valued_cells_to_moc(max_depth, uniq, values, 0_f64, 1_f64);
        let expect_nested_ranges = HpxRanges::new_unchecked(
            vec![0..12 << (2*29)]
        );

        assert_eq!(nested_ranges, expect_nested_ranges);
    }
}
