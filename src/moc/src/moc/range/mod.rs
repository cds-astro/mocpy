
use std::slice;
use std::ops::Range;
use std::vec::IntoIter;
use std::marker::PhantomData;
use std::convert::{TryInto, TryFrom};
use std::num::TryFromIntError;

use healpix::nested::{
  cone_coverage_approx_custom,
  elliptical_cone_coverage_custom,
  polygon_coverage,
  custom_polygon_coverage,
  zone_coverage,
  box_coverage,
  append_external_edge,
  bmoc::BMOC,
};
use healpix::sph_geom::ContainsSouthPoleMethod;

use crate::idx::Idx;
use crate::qty::{MocQty, Hpx, Time, Bounded};
use crate::elem::cell::Cell;
use crate::elemset::range::MocRanges;
use crate::moc::{
  HasMaxDepth, ZSorted, NonOverlapping, MOCProperties,
  RangeMOCIterator, RangeMOCIntoIterator
};
use crate::ranges::SNORanges;
use crate::moc::builder::{
  fixed_depth::{
    FixedDepthMocBuilder,
    OwnedOrderedFixedDepthCellsToRangesFromU64
  },
  maxdepth_range::RangeMocBuilder
};
use crate::moc::range::op::or::{or, OrRangeIter};
use crate::moc::range::op::minus::{minus, MinusRangeIter};
use crate::moc::range::op::and::{and, AndRangeIter};
use crate::moc::range::op::xor::xor;

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
  pub fn depth_max(&self) -> u8 {
    self.depth_max
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

  pub fn range_sum(&self) -> T {
    let mut sum = T::zero();
    for Range { start, end } in self.ranges.0.iter() {
      sum += *end - *start;
    }
    sum
  }

  pub fn coverage_percentage(&self) -> f64 {
    let rsum = self.range_sum();
    let tot = Q::upper_bound_exclusive();
    if T::N_BITS > 52 { // 52 = n mantissa bits in a f64
      // Divide by the same power of 2, dropping the LSBs
      let shift = (T::N_BITS - 52) as u32;
      rsum.unsigned_shr(shift);
      tot.unsigned_shr(shift);
    }
    rsum.cast_to_f64() / tot.cast_to_f64()
  }


  /// <=> from HEALPix map, i.e. from a list of HEALPic cell indices at the same depth
  pub fn from_fixed_depth_cells<I: Iterator<Item=T>>(
    depth: u8,
    cells_it: I,
    buf_capacity: Option<usize>
  ) -> Self {
    let mut builder = FixedDepthMocBuilder::new(depth, buf_capacity);
    for cell in cells_it {
      builder.push(cell)
    }
    builder.into_moc()
  }

  /// The value must be at the MOC depth
  pub fn contains_depth_max_val(&self, x: &T) -> bool {
    self.contains_val(&x.unsigned_shl(Q::shift_from_depth_max(self.depth_max) as u32))
  }

  pub fn contains_val(&self, x: &T) -> bool {
    self.ranges.contains_val(x)
  }

  pub fn append_fixed_depth_cells<I: Iterator<Item=T>>(
    self,
    depth: u8,
    cells_it: I,
    buf_capacity: Option<usize>
  ) -> Self {
    assert_eq!(depth, self.depth_max);
    let mut builder = FixedDepthMocBuilder::from(buf_capacity, self);
    for cell in cells_it {
      builder.push(cell)
    }
    builder.into_moc()
  }

  pub fn and(&self, rhs: &RangeMOC<T, Q>) -> RangeMOC<T, Q> {
    let depth_max = self.depth_max.max(rhs.depth_max);
    let ranges = self.ranges.intersection(&rhs.ranges);
    RangeMOC::new(depth_max, ranges)
  }
  pub fn intersection(&self, rhs: &RangeMOC<T, Q>) -> RangeMOC<T, Q> {
    self.and(rhs)
  }

  pub fn or(&self, rhs: &RangeMOC<T, Q>) -> RangeMOC<T, Q> {
    let depth_max = self.depth_max.max(rhs.depth_max);
    let ranges = self.ranges.union(&rhs.ranges);
    RangeMOC::new(depth_max, ranges)
  }
  pub fn union(&self, rhs: &RangeMOC<T, Q>) -> RangeMOC<T, Q> {
    self.or(rhs)
  }

  pub fn not(&self) -> RangeMOC<T, Q> {
    self.complement()
  }
  pub fn complement(&self) -> RangeMOC<T, Q> {
    RangeMOC::new(self.depth_max, self.ranges.complement())
  }

  pub fn xor(&self, rhs: &RangeMOC<T, Q>) -> RangeMOC<T, Q> {
    let depth_max = self.depth_max.max(rhs.depth_max);
    let ranges = xor((&self).into_range_moc_iter(), (&rhs).into_range_moc_iter()).collect();
    RangeMOC::new(depth_max, ranges)
  }

  pub fn minus(&self, rhs: &RangeMOC<T, Q>) -> RangeMOC<T, Q> {
    let depth_max = self.depth_max.max(rhs.depth_max);
    let ranges = minus((&self).into_range_moc_iter(), (&rhs).into_range_moc_iter()).collect();
    RangeMOC::new(depth_max, ranges)
  }

  pub fn degraded(&self, new_depth: u8) -> RangeMOC<T, Q> {
    let depth_max = self.depth_max.min(new_depth);
    let ranges = self.ranges.degraded(new_depth);
    RangeMOC::new(depth_max, ranges)
  }


  // CONTAINS: union that stops at first elem found
  // OVERLAP (=!CONTAINS on the COMPLEMENT ;) )

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
    /*println!("BMOC depth max: {}", bmoc.get_depth_max());
    for raw_cell in bmoc.iter() {
      let bmoc::Cell { raw_value, depth, hash, is_full } = bmoc.from_raw_value(*raw_cell);
      println!("depth: {}, idx: {}", depth, hash);
    }*/
    let shift = Hpx::<u64>::shift_from_depth_max(bmoc.get_depth_max());
    let mut ranges = bmoc.to_ranges();
    for range in ranges.iter_mut() {
      range.start <<= shift;
      range.end <<= shift;
    }
    // TODO: add a debug_assert! checking that the result is sorted!
    RangeMOC::new(bmoc.get_depth_max(), MocRanges::new_unchecked(ranges.to_vec()))
  }
}


impl<T: Idx> RangeMOC<T, Hpx<T>> {
  
  /// Add the MOC external border of depth `self.depth_max`.
  pub fn expanded(&self) -> Self {
    self.expanded_iter().into_range_moc()
  }

  pub fn expanded_iter(&self) -> OrRangeIter<
    T, 
    Hpx<T>,
    RangeRefMocIter<'_, T, Hpx<T>>, 
    OwnedOrderedFixedDepthCellsToRangesFromU64<T, Hpx<T>, IntoIter<u64>>
  > {
    let mut ext: Vec<u64> = Vec::with_capacity(10 * self.ranges.ranges().0.len()); // constant to be adjusted
    for Cell { depth, idx } in (&self).into_range_moc_iter().cells() {
      append_external_edge(depth, idx.to_u64(), self.depth_max - depth, &mut ext);
    }
    ext.sort_unstable(); // parallelize with rayon? It is the slowest part!!
    let ext_range_iter = OwnedOrderedFixedDepthCellsToRangesFromU64::new(self.depth_max, ext.into_iter());
    or((&self).into_range_moc_iter(), ext_range_iter)
  }

  pub fn contracted(&self) -> Self {
    self.contracted_iter().into_range_moc()
  }
  
  pub fn contracted_iter(&self) -> MinusRangeIter<
    T,
    Hpx<T>,
    RangeRefMocIter<'_, T, Hpx<T>>,
    OwnedOrderedFixedDepthCellsToRangesFromU64<T, Hpx<T>, IntoIter<u64>>
  > {
    let mut ext_of_complement: Vec<u64> = Vec::with_capacity(10 * self.ranges.ranges().0.len()); // constant to be adjusted
    for Cell { depth, idx } in (&self).into_range_moc_iter().not().cells() {
      append_external_edge(depth, idx.to_u64(), self.depth_max - depth, &mut ext_of_complement);
    }
    // _by(|a, b| a.flat_cmp::<Hpx<u64>>(&b))
    ext_of_complement.sort_unstable(); // parallelize with rayon? It is the slowest part!!
    let ext_range_iter = OwnedOrderedFixedDepthCellsToRangesFromU64::new(self.depth_max, ext_of_complement.into_iter());
    minus((&self).into_range_moc_iter(), ext_range_iter)
  }
  
  /// Returns this MOC external border
  pub fn external_border(&self) -> Self {
    /*minus(
      self.expanded_iter(),
      (&self).into_range_moc_iter()
    ).into_range_moc()*/
    self.external_border_iter().into_range_moc()
  }

  /// Returns this MOC external border
  pub fn external_border_iter(&self) -> MinusRangeIter<
    T,
    Hpx<T>,
    OrRangeIter<T, Hpx<T>, RangeRefMocIter<'_, T, Hpx<T>>, OwnedOrderedFixedDepthCellsToRangesFromU64<T, Hpx<T>, IntoIter<u64>>>,
    RangeRefMocIter<'_, T, Hpx<T>>
  > {
    minus(
      self.expanded_iter(),
      (&self).into_range_moc_iter()
    )
  }

  /// Returns this MOC internal border
  pub fn internal_border(&self) -> Self {
    let not = self.not();
    and(not.expanded_iter(), (&self).into_range_moc_iter()).into_range_moc()
  }

  pub fn internal_border_iter(&self) -> AndRangeIter<
    T, Hpx<T>,
    RangeMocIter<T, Hpx<T>>,
    RangeRefMocIter<'_, T, Hpx<T>>
  > {
    let left = self.not().expanded();
    and(left.into_range_moc_iter(), (&self).into_range_moc_iter())
  }
  
}


fn from<T: Idx + TryFrom<u64, Error=TryFromIntError>>(range_moc: RangeMOC<u64, Hpx<u64>>) -> RangeMOC<T, Hpx<T>> {
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
  
  /// # Panics
  ///   If a `lat` is **not in** `[-pi/2, pi/2]`, this method panics.
  pub fn from_coos<I: Iterator<Item=(f64, f64)>>(depth: u8, coo_it: I, buf_capacity: Option<usize>) -> Self {
    let layer = healpix::nested::get(depth);
    Self::from_fixed_depth_cells(
      depth,
      coo_it.map(move |(lon_rad, lat_rad)| layer.hash(lon_rad, lat_rad)),
      buf_capacity
    )
  }

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
  /// - `lon` the longitude of the center of the box, in radians
  /// - `lat` the latitude of the center of the box, in radians
  /// - `a` the semi-major axis of the box (half the box width), in radians
  /// - `b` the semi-minor axis of the box (half the box height), in radians
  /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
  ///
  /// # Panics
  /// - if `a` not in `]0, pi/2]`
  /// - if `b` not in `]0, a]`
  /// - if `pa` not in `[0, pi[`
  /// 
  pub fn from_box(lon: f64, lat: f64, a: f64, b: f64, pa: f64, depth: u8) -> Self {
    Self::from(box_coverage(depth, lon, lat, a, b, pa))
  }
  
  /// # Input
  /// - `vertices` the list of vertices (in a slice) coordinates, in radians
  ///              `[(lon, lat), (lon, lat), ..., (lon, lat)]`
  /// - `depth`: the MOC depth
  pub fn from_polygon(vertices: &[(f64, f64)], complement: bool, depth: u8) -> Self {
    Self::from(
      if !complement {
        polygon_coverage(depth, vertices, true)
      } else {
        custom_polygon_coverage(depth, vertices, &ContainsSouthPoleMethod::DefaultComplement, true)
      }
    )
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
  
  // BORDER = NOT(SELF)+EXPAND && SELF

  /* Perform UNIONS
  pub fn from_fixed_radius_cones
  pub fn from_multi_cones
  pub fn from_multi_elliptical_cones*/

  /*
  /// Add the MOC external border of depth `self.depth_max`.
  pub fn expanded(&self) -> Self {
    self.expanded_iter().into_range_moc()
  }

  pub fn expanded_iter(&self) -> OrRangeIter<u64, Hpx<u64>,
    RangeRefMocIter<'_, u64, Hpx<u64>>, OwnedOrderedFixedDepthCellsToRanges<u64, Hpx<u64>>> {
    let mut ext: Vec<u64> = Vec::with_capacity(10 * self.ranges.ranges().0.len()); // constant to be adjusted
    for Cell { depth, idx } in (&self).into_range_moc_iter().cells() {
      append_external_edge(depth, idx, self.depth_max - depth, &mut ext);
    }
    ext.sort_unstable(); // parallelize with rayon? It is the slowest part!!
    let ext_range_iter = OwnedOrderedFixedDepthCellsToRanges::new(self.depth_max, ext.into_iter());
    or((&self).into_range_moc_iter(), ext_range_iter)
  }

  /// Returns this MOC external border
  pub fn external_border(&self) -> Self {
    self.external_border_iter().into_range_moc()
  }

  /// Returns this MOC external border
  pub fn external_border_iter(&self) -> MinusRangeIter<
    u64,
    Hpx<u64>,
    OrRangeIter<u64, Hpx<u64>, RangeRefMocIter<'_, u64, Hpx<u64>>, OwnedOrderedFixedDepthCellsToRanges<u64, Hpx<u64>>>,
    RangeRefMocIter<'_, u64, Hpx<u64>>
  > {
    minus(
      self.expanded_iter(),
      (&self).into_range_moc_iter()
    )
  }

  /// Returns this MOC internal border
  pub fn internal_border(&self) -> Self {
    let not = self.not();
    and(not.expanded_iter(), (&self).into_range_moc_iter()).into_range_moc()
  }

  pub fn internal_border_iter(&self) -> AndRangeIter<
    u64, Hpx<u64>,
    RangeMocIter<u64, Hpx<u64>>,
    RangeRefMocIter<'_, u64, Hpx<u64>>
  > {
    let left = self.not().expanded();
    and(left.into_range_moc_iter(), (&self).into_range_moc_iter())
  }*/

}


impl<T: Idx> RangeMOC<T, Time<T>> {

  pub fn from_microsec_since_jd0<I>(depth: u8, it: I, buf_capacity: Option<usize>) -> Self
    where
      I: Iterator<Item=u64>
  {
    let shift = Time::<T>::shift_from_depth_max(depth) as u32;
    let mut builder = FixedDepthMocBuilder::new(depth, buf_capacity);
    for t in it {
      builder.push(T::from_u64_idx(t).unsigned_shr(shift));
    }
    builder.into_moc()
  }

  pub fn from_microsec_ranges_since_jd0<I>(depth: u8, it: I, buf_capacity: Option<usize>) -> Self
    where
      I: Iterator<Item=Range<u64>>
  {
    let mut builder = RangeMocBuilder::new(depth, buf_capacity);
    for range in it {
      builder.push(T::from_u64_idx(range.start)..T::from_u64_idx(range.end));
    }
    builder.into_moc()
  }

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
