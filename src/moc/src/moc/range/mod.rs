
use std::slice;
use std::ops::Range;
use std::vec::{IntoIter};
use std::marker::PhantomData;

use healpix::nested::{
  cone_coverage_approx_custom,
  elliptical_cone_coverage_custom,
  polygon_coverage,
  zone_coverage,
  bmoc::BMOC,
};

use crate::idx::Idx;
use crate::qty::{MocQty, Hpx};
use crate::elemset::range::MocRanges;
use crate::moc::{
  HasMaxDepth, ZSorted, NonOverlapping, MOCProperties,
  RangeMOCIterator, RangeMOCIntoIterator
};
use std::convert::{TryInto, TryFrom};
use std::num::TryFromIntError;

pub mod op;

/// A MOC made of (ordered and non-overlaping) ranges.
#[derive(Debug, Clone)]
pub struct RangeMOC<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  ranges: MocRanges<T, Q>
}
impl<T: Idx, Q: MocQty<T>> RangeMOC<T, Q> {
  pub fn new(depth_max: u8, ranges: MocRanges<T, Q>) -> Self {
    Self {depth_max, ranges }
  }
  /// Returns the number of ranges the MOC contains
  pub fn len(&self) -> usize {
    self.ranges.0.0.len()
  }
  pub fn is_empty(&self) -> bool { self.len() == 0 }
  pub fn moc_ranges(&self) -> &MocRanges<T, Q> {
    &self.ranges
  }
  pub fn into_moc_ranges(self) -> MocRanges<T, Q> {
    self.ranges
  }

  // pub fn and(&self, RangeMocIntoIterator) !!

  // pub fn owned_and() -> RangeMOC<T, Q> { }
  // pub fn lazy_and() -> 
  
  
  /*pub fn into_range_moc_iter(self) -> LazyRangeMOCIter<T, Q, IntoIter<Range<T>>> {
    LazyRangeMOCIter::new(self.depth_max, self.ranges.0.0.into_iter())
  }*/

  /*pub fn range_moc_iter(&self) -> LazyRangeMOCVecIter<'_, H> {
    LazyRangeMOCVecIter::new(self.depth_max, self.ranges.iter())
  }*/
  /*pub fn into_cells_iter(self) -> CellMOCIteratorFromRanges<T, Q, Self> {
    CellMOCIteratorFromRanges::new(self)
  }*/
  /*pub fn to_cells_iter(&self) -> CellMOCIteratorFromRanges<T, Q, Self> {
    CellMOCIteratorFromRanges::new(self)
  }*/
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for RangeMOC<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for RangeMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for RangeMOC<T, Q> { }

impl From<BMOC> for RangeMOC<u64, Hpx<u64>> {
  fn from(bmoc: BMOC) -> Self {
    let ranges = bmoc.to_ranges();
    // TODO: add a debug_assert! checking that the result is sorted!
    RangeMOC::new(bmoc.get_depth_max(), MocRanges::new_unchecked(ranges.to_vec()))
  }
}


fn from<T: Idx + TryFrom<u64,Error=TryFromIntError>>(range_moc: RangeMOC<u64, Hpx<u64>>) -> RangeMOC<T, Hpx<T>> {
  let depth_max= range_moc.depth_max;
  let ranges = range_moc.ranges.0;
  let shift = u64::N_BITS - T::N_BITS;
  let ranges: Vec<Range<T>> = ranges.0.iter()
      .map(|Range { start, end}| (start >> shift).try_into().unwrap()..(end >> shift).try_into().unwrap())
      .collect();
  RangeMOC::new(depth_max, MocRanges::new_unchecked(ranges))
}

impl From<RangeMOC<u64, Hpx<u64>>> for RangeMOC<u32, Hpx<u32>> {
  fn from(range_moc: RangeMOC<u64, Hpx<u64>>) -> Self {
    assert!(range_moc.depth_max < 14);
    from(range_moc)
  }
}

impl From<RangeMOC<u64, Hpx<u64>>> for RangeMOC<u16, Hpx<u16>> {
  fn from(range_moc: RangeMOC<u64, Hpx<u64>>) -> Self {
    assert!(range_moc.depth_max < 6);
    from(range_moc)
  }
}

impl RangeMOC<u64, Hpx<u64>> {

  /// # Input
  /// - `cone_lon` the longitude of the center of the cone, in radians
  /// - `cone_lat` the latitude of the center of the cone, in radians
  /// - `cone_radius` the radius of the cone, in radians
  /// - `depth`: the MOC depth
  /// - `delta_depth` the difference between the MOC depth and the depth at which the computations
  ///   are made (should remain quite small).
  ///
  /// # Panics
  /// If this layer depth + `delta_depth` > the max depth (i.e. 29)
  pub fn from_cone(lon: f64, lat: f64, radius: f64, depth: u8, delta_depth: u8) -> Self {
    Self::from(cone_coverage_approx_custom(depth, delta_depth, lon, lat, radius))
  }

  /// # Input
  /// - `lon` the longitude of the center of the elliptical cone, in radians
  /// - `lat` the latitude of the center of the elliptical cone, in radians
  /// - `a` the semi-major axis of the elliptical cone, in radians
  /// - `b` the semi-minor axis of the elliptical cone, in radians
  /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
  /// - `depth`: the MOC depth
  /// - `delta_depth` the difference between the MOC depth and the depth at which the computations
  ///   are made (should remain quite small).
  ///
  /// # Panics
  /// - if the semi-major axis is > PI/2
  /// - if this layer depth + `delta_depth` > the max depth (i.e. 29)
  pub fn from_elliptical_cone(lon: f64, lat: f64, a: f64, b: f64, pa: f64, depth: u8, delta_depth: u8) -> Self {
    Self::from(elliptical_cone_coverage_custom(depth, delta_depth, lon, lat, a, b, pa))
  }

  /// # Input
  /// - `vertices` the list of vertices (in a slice) coordinates, in radians
  ///              `[(lon, lat), (lon, lat), ..., (lon, lat)]`
  /// - `depth`: the MOC depth
  pub fn from_polygon(vertices: &[(f64, f64)], depth: u8) -> Self {
    Self::from(polygon_coverage(depth, vertices, true))
  }

  /// # Input
  /// - `lon_min` the longitude of the bottom left corner
  /// - `lat_min` the latitude of the bottom left corner
  /// - `lon_max` the longitude of the upper left corner
  /// - `lat_max` the latitude of the upper left corner
  /// - `depth`: the MOC depth
  ///
  /// # Remark
  /// - If `lon_min > lon_max` then we consider that the zone crosses the primary meridian.
  /// - The north pole is included only if `lon_min == 0 && lat_max == pi/2`
  ///
  /// # Panics
  /// * if `lon_min` or `lon_max` not in `[0, 2\pi[`
  /// * if `lat_min` or `lat_max` not in `[-\pi/2, \pi/2[`
  /// * `lat_min >= lat_max`.
  pub fn from_zone(lon_min: f64, lat_min: f64, lon_max: f64, lat_max: f64, depth: u8) -> Self {
    Self::from(zone_coverage(depth, lon_min, lat_min, lon_max, lat_max))
  }

  /*
  expand
  c.f. healpix: moc/mod.rs L455
  pub fn expand(self) -> OrMocIter<u64, MOCIteratorFromRanges<u64,LazyRangeMoc>, MergeMoc> {
    let mut ext: Vec<u64> = Vec::with_capacity(10 * self.ranges.len()); // constant to be adjusted
    for HpxCell { depth, hash } in MOCIteratorFromRanges::new(self.range_moc_iter()) {
      append_external_edge(depth, hash, self.depth_max - depth, &mut ext);
    }
    ext.sort_unstable(); // parallelize with rayon? It is the slowest part!!
    ext.dedup();
    let depth = self.depth_max;
    MOCIteratorFromRanges::new(self.into_range_moc_iter())
        .or(MergeIter::new(FlatLazyMOCIter::new(depth, ext.into_iter())
        )
        )
  }
  */

  /*
  pub fn from_cells => CellMOC -> Ranges
  pub fn from_coos  => CellMOC -> Ranges
  */

  /* Perform UNIONS
  pub fn from_fixed_radius_cones
  pub fn from_multi_cones
  pub fn from_multi_elliptical_cones*/

}


/// Iterator taking the ownership of the `RangeMOC` it iterates over.
pub struct RangeMocIter<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: IntoIter<Range<T>>,
  last: Option<Range<T>>,
  _qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for RangeMocIter<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> Iterator for RangeMocIter<T, Q> {
  type Item = Range<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
  // Declaring size_hint, a 'collect' can directly allocate the right number of elements
  fn size_hint(&self) -> (usize, Option<usize>) {
    self.iter.size_hint()
  }
}
impl<T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for RangeMocIter<T, Q> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    self.last.as_ref()
  }
}
impl<T: Idx, Q: MocQty<T>> RangeMOCIntoIterator<T> for RangeMOC<T, Q> {
  type Qty = Q;
  type IntoRangeMOCIter = RangeMocIter<T, Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter {
    let l = self.ranges.0.0.len();
    let last: Option<Range<T>> = if l > 0 {
      Some(self.ranges.0.0[l - 1].clone())
    } else {
      None
    };
    RangeMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.into_vec().into_iter(),
      last,
      _qty: PhantomData
    }
  }
}

/// Iterator borrowing the `RangeMOC` it iterates over.
pub struct RangeRefMocIter<'a, T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: slice::Iter<'a, Range<T>>,
  last: Option<Range<T>>,
  _qty: PhantomData<Q>,
}
impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for RangeRefMocIter<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<'a, T: Idx, Q: MocQty<T>> ZSorted for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> Iterator for RangeRefMocIter<'a, T, Q> {
  type Item = Range<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().cloned()
  }
  // Declaring size_hint, a 'collect' can directly allocate the right number of elements
  fn size_hint(&self) -> (usize, Option<usize>) {
    self.iter.size_hint()
  }
}
impl<'a, T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for RangeRefMocIter<'a, T, Q> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    self.last.as_ref()
  }
}
impl<'a, T: Idx, Q: MocQty<T>> RangeMOCIntoIterator<T> for &'a RangeMOC<T, Q> {
  type Qty = Q;
  type IntoRangeMOCIter = RangeRefMocIter<'a, T, Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter {
    let l = self.ranges.0.0.len();
    let last: Option<Range<T>> = if l > 0 {
      Some(self.ranges.0.0[l - 1].clone())
    } else {
      None
    };
    RangeRefMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.iter(),
      last,
      _qty: PhantomData
    }
  }
}
