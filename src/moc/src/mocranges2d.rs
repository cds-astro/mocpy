
use std::ops::Range;
use std::marker::PhantomData;

use rayon::prelude::*;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::elemset::range::MocRanges;
use crate::ranges::{
  Ranges,
  ranges2d::{Ranges2D, SNORanges2D}
};

// Declaration of the ST-MOC type mage in hpxranges2d
// pub type TimeSpaceMoc<T, S> = Moc2DRanges::<T, Time<T>, S, Hpx<S>>;

// only removing the depth() from Ranges2D<T, S> and let it in Moc2DRanges<TT, T, ST, S>
#[derive(Debug)]
pub struct Moc2DRanges<TT, T, ST, S>
  where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
    S: MocQty<ST>,
{
  pub ranges2d: Ranges2D<TT, ST>,
  _t_type: PhantomData<T>,
  _s_type: PhantomData<S>,
}

impl<'a, TT, T, ST, S> SNORanges2D<'a, TT, ST> for Moc2DRanges<TT, T, ST, S>
  where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
    S: MocQty<ST>,
{

  fn make_consistent(mut self) -> Self {
    self.ranges2d = self.ranges2d.make_consistent();
    self
  }

  fn is_empty(&self) -> bool {
    self.ranges2d.is_empty()
  }

  fn contains(&self, time: TT, range: &Range<ST>) -> bool {
    self.ranges2d.contains(time, range)
  }

  fn union(&self, other: &Self) -> Self {
    self.ranges2d.union(&other.ranges2d).into()
  }

  fn intersection(&self, other: &Self) -> Self {
    self.ranges2d.intersection(&other.ranges2d).into()
  }

  fn difference(&self, other: &Self) -> Self {
    self.ranges2d.difference(&other.ranges2d).into()
  }
}

impl <TT, T, ST, S> From<Ranges2D<TT, ST>> for Moc2DRanges<TT, T, ST, S>
  where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
    S: MocQty<ST> {
  fn from(ranges2d: Ranges2D<TT, ST>) -> Self {
    Moc2DRanges {
      ranges2d,
      _t_type: PhantomData,
      _s_type: PhantomData
    }
  }
}

/*impl <TT, T, ST, S> From<Moc2DRanges<TT, T, ST, S>> for Ranges2D<TT, ST>
  where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
    S: MocQty<ST> {
  fn from(mocranges2d: Moc2DRanges<TT, T, ST, S>) -> Self {
    mocranges2d.ranges2d
  }
}*/


impl<TT, T, ST, S> Moc2DRanges<TT, T, ST, S>
  where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
    S: MocQty<ST>,
{
  /// Creates a 2D coverage
  ///
  /// # Arguments
  ///
  /// * `t` - A set of ranges constituing the first dimension. This stores
  ///   usually quantities such as times, redshifts or proper motions.
  /// * `s` - A set of 1D coverage constituing the second dimension. This stores
  ///   usually space informations such as HEALPix cell indices under the nested format.
  ///
  /// # Precondition
  ///
  /// ``t`` and ``s`` must have the same length.
  pub fn new(t: Vec<Range<TT>>, s: Vec<Ranges<ST>>) -> Moc2DRanges<TT, T, ST, S> {
    Moc2DRanges {
      ranges2d: Ranges2D::new(t, s),
      _t_type: PhantomData,
      _s_type: PhantomData
    }
  }

  /// Compute the smallest possible depth of the coverage
  ///
  /// # Returns
  ///
  /// A tuple containing two values:
  ///
  /// * The maximum depth along the `T` axis
  /// * The maximum depth along the `S` axis
  ///
  /// # Info
  ///
  /// If the `NestedRanges2D<T, S>` is empty, the depth returned
  /// is set to (0, 0)
  pub fn compute_min_depth(&self) -> (u8, u8) {
    let y = self.ranges2d.y
      .par_iter()
      // Compute the depths of the Ranges<S>
      .map(|ranges| MocRanges::<ST, S>::compute_min_depth_gen(ranges))
      // Get the max of these depths
      .max()
      // If there are no ranges, the max depth
      // along the second dimension is set to 0
      .unwrap_or(0);

    // The computation is very light (logical OR), so I wonder about the the cost (overhead)
    // of the parallelization here (except for very large MOCs)...
    let x = T::compute_min_depth(
      self.ranges2d.x.par_iter()
        // Perform a logical 'or' between (upper and lower bounds of) all indices of the first dimension
        // then look at the trailing zeros (in the compute_min_depth method)
        .fold_with(TT::zero(), |acc, range| acc | range.start | range.end)
        .reduce(TT::zero, |a, b| a | b)
    );

    (x, y)
  }

}

impl<TT, T, ST, S> PartialEq for Moc2DRanges<TT, T, ST, S>
  where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
    S: MocQty<ST>,
{
  fn eq(&self, other: &Self) -> bool {
    self.ranges2d.eq(&other.ranges2d)
  }
}

#[cfg(test)]
mod tests {
  use std::ops::Range;

  use crate::idx::Idx;
  use crate::ranges::Ranges;
  use crate::ranges::ranges2d::{SNORanges2D};
  use crate::qty::{Time, Hpx};
  use crate::mocranges2d::Moc2DRanges;
  use crate::hpxranges2d::{TimeSpaceMoc, HpxRanges2D};

  // Tests are the same as (a sub-part of) range2d test.
  // So  we basically test with the decorator (algo are already tested in rande2d).

  type TimeSpaceRanges<T, S> = Moc2DRanges::<T, Time<T>, S, Hpx<S>>;

  fn time_space_moc<T: Idx, S: Idx>(ranges: Moc2DRanges::<T, Time<T>, S, Hpx<S>>) -> TimeSpaceMoc<T, S> {
    HpxRanges2D(ranges)
  }

  fn new_time_space_moc<T: Idx, S: Idx>(t: Vec<Range<T>>, s: Vec<Ranges<S>>) -> TimeSpaceMoc<T, S> {
    time_space_moc(TimeSpaceRanges::<T, S>::new(t, s).make_consistent())
  }

  #[test]
  fn merge_overlapping_ranges() {
    let t: Vec<Range<u64>> = vec![0..15, 0..15, 15..30, 30..45, 15..30];
    let s = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
    ];
    let coverage = new_time_space_moc(t, s);

    let t_expect = vec![0..15, 15..45];
    let s_expect = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
    ];
    let coverage_expect = new_time_space_moc(t_expect, s_expect);
    assert_eq!(coverage, coverage_expect);
  }

  // Overlapping time ranges configuration:
  // xxxxxxxxxxx
  // xxxx-------
  #[test]
  fn remove_different_length_time_ranges() {
    let t: Vec<Range<u64>> = vec![0..7, 0..30];
    let s = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
    ];
    let coverage = new_time_space_moc(t, s);

    let t_expect = vec![0..7, 7..30];
    let s_expect = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
      Ranges::<u64>::new_unchecked(vec![16..21])
    ];
    let coverage_expect = new_time_space_moc(t_expect, s_expect);
    assert_eq!(coverage, coverage_expect);
  }

  // Overlapping time ranges configuration:
  // xxxxxxxxxxx
  // ----xxxx---
  #[test]
  fn remove_different_length_time_ranges2() {
    let t: Vec<Range<u64>> = vec![0..30, 2..10];
    let s = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
    ];
    let coverage = new_time_space_moc(t, s);

    let t_expect = vec![0..2, 2..10, 10..30];
    let s_expect = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18])
    ];
    let coverage_expect = new_time_space_moc(t_expect, s_expect);
    assert_eq!(coverage, coverage_expect);
  }

  // Overlapping time ranges configuration:
  // xxxxxxx----
  // ----xxxxxxx
  #[test]
  fn remove_different_length_time_ranges3() {
    let t: Vec<Range<u64>> = vec![0..5, 2..10];
    let s = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
    ];
    let coverage = new_time_space_moc(t, s);

    let t_expect = vec![0..2, 2..5, 5..10];
    let s_expect = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
      Ranges::<u64>::new_unchecked(vec![16..21])
    ];
    let coverage_expect = new_time_space_moc(t_expect, s_expect);
    assert_eq!(coverage, coverage_expect);
  }

  // Overlapping time ranges configuration:
  // xxxxxxxxxxx
  // ----xxxxxxx
  #[test]
  fn remove_different_length_time_ranges4() {
    let t: Vec<Range<u64>> = vec![0..30, 10..30];
    let s = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
    ];
    let coverage = new_time_space_moc(t, s);

    let t_expect = vec![0..10, 10..30];
    let s_expect = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
    ];
    let coverage_expect = new_time_space_moc(t_expect, s_expect);
    assert_eq!(coverage, coverage_expect);
  }
  // No overlapping time ranges
  // xxxxxx----
  // ------xxxx
  #[test]
  fn remove_different_length_time_ranges5() {
    let t: Vec<Range<u64>> = vec![0..5, 5..20];
    let s = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
    ];
    let coverage = new_time_space_moc(t, s);

    let t_expect = vec![0..5, 5..20];
    let s_expect = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![16..21])
    ];
    let coverage_expect = new_time_space_moc(t_expect, s_expect);
    assert_eq!(coverage, coverage_expect);
  }

  #[test]
  fn merge_overlapping_ranges_2() {
    let t: Vec<Range<u64>> = vec![0..15, 0..15, 15..30, 30..45, 15..30];
    let s = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![0..4]),
      Ranges::<u64>::new_unchecked(vec![16..21]),
      Ranges::<u64>::new_unchecked(vec![16..21, 25..26]),
    ];
    let coverage = new_time_space_moc(t, s);

    let t_expect = vec![0..15, 15..30, 30..45];
    let s_expect = vec![
      Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
      Ranges::<u64>::new_unchecked(vec![0..4, 16..21, 25..26]),
      Ranges::<u64>::new_unchecked(vec![16..21])
    ];
    let coverage_expect = new_time_space_moc(t_expect, s_expect);
    assert_eq!(coverage, coverage_expect);
  }

}