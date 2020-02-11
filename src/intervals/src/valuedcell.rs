use num::{Num, Integer, PrimInt, One};
use ndarray::Array1;
use rayon::prelude::*;

use std::cmp::Ordering::Equal;
use std::ops::Range;

use crate::bounded::{Bounded, NestedRange};
use crate::nestedranges::NestedRanges;

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
    max_depth: u32,
    uniq: Array1<T>,
    values: Array1<V>,
    cumul_from: V,
    cumul_to: V
    ) -> NestedRanges<T>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      V: Num + PartialOrd + DivBy<T, Output=V> + Copy + Send + Sync + std::fmt::Debug {
    let mut valued_uniq_sorted: Vec<(T, V, V)> = uniq.iter().zip(values.iter())
        .map(|(uniq, val)| {
            let (depth, _icell) = T::pix_depth(*uniq);
            let n_sub_cells = T::one().unsigned_shl((max_depth - depth) << 1);
            (*uniq, *val, val.div_by(n_sub_cells))
        })
    .collect::<Vec<(T, V, V)>>();
    // We use b.comp(a) instead of a.cmp(b) to get the DESC order
    valued_uniq_sorted.par_sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Equal));
    let mut result: Vec<Range<T>> = Vec::with_capacity(valued_uniq_sorted.len());

    let mut i = 0_usize;
    let mut acc = V::zero();
    while i < valued_uniq_sorted.len() && acc.add(valued_uniq_sorted[i].1) <= cumul_from {
        acc = acc.add(valued_uniq_sorted[i].1);
        result.push(T::uniq_to_range(valued_uniq_sorted[i].0));
        i += 1;
    }
    if i < valued_uniq_sorted.len() && acc < cumul_from {
        let (depth, icell) = T::pix_depth(valued_uniq_sorted[i].0);
        result = recursive_descent(
            depth, icell, max_depth,
            valued_uniq_sorted[i].1, cumul_from.sub(acc), result);
        i += 1;
    }
    while i < valued_uniq_sorted.len() && acc.add(valued_uniq_sorted[i].1) <= cumul_to {
        acc = acc.add(valued_uniq_sorted[i].1);
        result.push(T::uniq_to_range(valued_uniq_sorted[i].0));
        i += 1;
    }
    if i < valued_uniq_sorted.len() && acc < cumul_to {
        let (depth, icell) = T::pix_depth(valued_uniq_sorted[i].0);
        result = recursive_descent(
            depth, icell, max_depth,
            valued_uniq_sorted[i].1, cumul_to.sub(acc), result);
    }
    NestedRanges::new(result).make_consistent()
}

fn recursive_descent<T, V>(
    depth: u32, ipix: T, max_depth: u32,
    cell_val: V, mut target_val: V,
    mut result: Vec<Range<T>>
    ) -> Vec<Range<T>>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      V: Num + PartialOrd + DivBy<T, Output=V> + Copy + Send + Sync {
    if depth == max_depth {
        if cell_val <= target_val {
            let rng: NestedRange<T> = (depth as u8, ipix).into();
            result.push(rng.0);
        }
    } else if target_val > V::zero() {
        let four = T::one().unsigned_shl(2);

        let subcell_val = cell_val.div_by(four);
        let depth = depth + 1;
        let ipix = ipix << 2;
        let mut i = T::zero();
        while i < four && target_val.sub(subcell_val) >= V::zero() {
            let rng: NestedRange<T> = (depth as u8, ipix + i).into();
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
    use crate::nestedranges::NestedRanges;
    use ndarray::prelude::*;
    use crate::valuedcell::valued_cells_to_moc;
    use crate::bounded::Bounded;
    use std::f64;
    use std::u64;
    use num::PrimInt;
    use crate::valuedcell::DivBy;

    #[test]
    fn test_single_uniq() {
        let uniq = array![4];
        let values = array![1_f64];

        let max_depth = 2;

        let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 0.25_f64);
        let tdd = (<u64>::MAXDEPTH as u32 - max_depth) << 1;
        let expect_nested_ranges = NestedRanges::new(
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

        let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 1_f64);
        let tdd = (<u64>::MAXDEPTH as u32 - max_depth) << 1;
        let expect_nested_ranges = NestedRanges::new(
            vec![]
        );

        assert_eq!(nested_ranges, expect_nested_ranges);
    }

    #[test]
    fn test_full_space() {
        let uniq = array![4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        let values = array![0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0_f64, 0_f64];

        let max_depth = 2;

        let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 1_f64);
        let tdd = (<u64>::MAXDEPTH as u32 - max_depth) << 1;
        let expect_nested_ranges = NestedRanges::new(
            vec![0..12 << (2*29)]
        );

        assert_eq!(nested_ranges, expect_nested_ranges);
    }
}
