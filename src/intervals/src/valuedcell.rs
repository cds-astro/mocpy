use num::{Num, Integer, PrimInt, One};
use ndarray::Array1;
use rayon::prelude::*;

use std::cmp::Ordering::Equal;
use std::ops::{Div, Range};

use crate::bounded::Bounded;
use crate::nestedranges::NestedRanges;

/// Creates a MOC from the given list of uniq cells numbers according to the value they contains.
/// We assume that the value is directly proportional to the covered area (like a flux or a probability).
/// Limits are put to select an area having a cumulative value ranging from a given lower limit
/// to a given upper limit.
/// An example is the selection of a region having between 10 and 90 percent of a flux, or
/// an 90 percent completeness.
/// 
/// # Args
/// * `max_depth`: the largest depth of the output MOC, which must be larger or equals to the largest
/// depth in the `uniq` values
/// * `uniq`: the list of uniq cells (i.e. values encoding both the HEALPix depth and cell number)
/// * `values`: values associated to each uniq. The table must have the same size as the `uniq` array.
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
        V: Num + PartialOrd + Div<T, Output=V> + Copy + Send + Sync {
  assert_eq!(uniq.len(), values.len());
  assert!(cumul_from < cumul_to);
  let mut valued_uniq_sorted: Vec<(T, V, V)> = uniq.iter().zip(values.iter())
    .map(|(uniq, val)| {
      let (depth, _icell) = T::pix_depth(*uniq);
      let n_sub_cells = T::one().unsigned_shl((max_depth - depth) << 1);
      (*uniq, *val, val.div(n_sub_cells))
    })
    .collect::<Vec<(T, V, V)>>();
  // We use b.comp(a) instead of a.cmp(b) to get the DESC order
  valued_uniq_sorted.par_sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Equal)); 
  let mut result: Vec<Range<T>> = Vec::with_capacity(valued_uniq_sorted.len());
  if valued_uniq_sorted.len() < 1 {
    return NestedRanges::new(result);
  }
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
  NestedRanges::new(result)
}

fn recursive_descent<T, V>(
  depth: u32, ipix: T, max_depth: u32,
  cell_val: V, mut target_val: V,
  mut result: Vec<Range<T>>
) -> Vec<Range<T>>
  where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
        V: Num + PartialOrd + Div<T, Output=V> + Copy + Send + Sync {
  let four = T::one().unsigned_shl(2);
  if depth == max_depth {
    if cell_val <= target_val {
      result.push(T::uniq_to_range(T::to_uniq(depth, ipix)));
    }
  } else if target_val > V::zero() {
    let subcell_val = cell_val.div(four);
    let depth = depth + 1;
    let ipix = ipix << 2;
    let mut i = T::zero();
    while i < four && target_val.sub(subcell_val) >= V::zero() {
      result.push(T::uniq_to_range(T::to_uniq(depth, ipix + i)));
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
