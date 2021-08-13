
use std::ops::Range;
use std::cmp::Ordering::Equal;

use num::{Num, One};

use crate::idx::Idx;
use crate::qty::Hpx;

use super::range::HpxRange;
use super::super::elemset::range::HpxRanges;

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
pub fn valued_cells_to_moc<'a, T, V, I1, I2>(
    max_depth: u8,
    uniq: I1,
    values: I2,
    cumul_from: V,
    cumul_to: V
) -> HpxRanges<T>
    where 
      T: Idx,
      V: 'static + Num + PartialOrd + DivBy<T, Output=V> + Copy + Send + Sync + std::fmt::Debug,
      I1: Iterator<Item=&'a T>,
      I2: Iterator<Item=&'a V>,
{
    let mut valued_uniq_sorted: Vec<(T, V, V)> = uniq.zip(values)
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
            i += One::one();
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

    use crate::qty::{MocQty, Hpx};
    use crate::elem::valuedcell::valued_cells_to_moc;
    use crate::elemset::range::HpxRanges;

    #[test]
    fn test_single_uniq() {
        let uniq = vec![4];
        let values = vec![1_f64];

        let max_depth = 2;

        // let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 0.25_f64);
        let nested_ranges = valued_cells_to_moc(max_depth, uniq.iter(), values.iter(), 0_f64, 0.25_f64);

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
        let uniq = vec![];
        let values = vec![];

        let max_depth = 2;

        // let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 1_f64);
        let nested_ranges = valued_cells_to_moc(max_depth, uniq.iter(), values.iter(), 0_f64, 1_f64);
        let expect_nested_ranges = HpxRanges::default();

        assert_eq!(nested_ranges, expect_nested_ranges);
    }

    #[test]
    fn test_full_space() {
        let uniq = vec![4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        let values = vec![0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0.1_f64, 0_f64, 0_f64];

        let max_depth = 2;

        // let nested_ranges = valued_cells_to_moc::<u64, f64>(max_depth, uniq, values, 0_f64, 1_f64);
        let nested_ranges = valued_cells_to_moc(max_depth, uniq.iter(), values.iter(), 0_f64, 1_f64);
        let expect_nested_ranges = HpxRanges::new_unchecked(
            vec![0..12 << (2*29)]
        );

        assert_eq!(nested_ranges, expect_nested_ranges);
    }
}
