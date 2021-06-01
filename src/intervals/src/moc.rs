
use std::slice;
use std::ops::Range;
use std::vec::{Vec, IntoIter};
use std::marker::PhantomData;

use crate::mocrange::MocRange;
use crate::ranges::Idx;
use crate::mocqty::MocQty;
use crate::mocranges::MocRanges;
use crate::mocell::{MocCell, Cell};


/// Returns the maximum depth of an item the implementor contains.
pub trait HasMaxDepth {
  fn depth_max(&self) -> u8;
}

/// Marker trait telling that something the implementor is sorted
pub trait Sorted {}

/// Marker trait telling that something the implementor contains non-overlapping items.
pub trait NonOverlapping {}

pub trait MOCProperties: HasMaxDepth + Sorted + NonOverlapping { }

/*pub trait MOC: MOCProperties {
  type Qty: MocQty<T>;
}*/

// HierarchMOC
// RangesMOC
// CellsMOC



/*pub trait CellMOCIterator<T, Q>: Sized + MOCProperties + Iterator<Item=MocCell<T, Q>>
  where
    T: Idx,
    Q: MocQty<T>,
{

}
impl<T, Q, R> CellMOCIterator<T, Q> for R
  where
    T: Idx,
    Q: MocQty<T>,
    R: Sized + MOCProperties + Iterator<Item=MocCell<T, Q>> {}
*/

pub trait CellMOCIterator<T: Idx>: Sized + MOCProperties + Iterator<Item=Cell<T>> {
  type Qty: MocQty<T>;
}

pub trait CellMOCIntoIterator<T: Idx>: Sized + MOCProperties + IntoIterator<Item=Cell<T>> {
  type Qty: MocQty<T>;
}

// RangeMOCIterator<T, Q, B>: Sized + MOCProperties + Iterator<Item=B>
// with B: Borrow<MocRange<T, Q>>
// Don't do this, we will clone as iterating if necessary

/*
/// Here we define a MOC as a simple sequence of ranges having the following properties:
/// * the range bounds are expressed at Q::MAX_DEPTH
/// * the list contains mutually exclusive ranges
/// * the list is sorted from the smaller to the largest bounds
pub trait RangeMOCIterator<T, Q>: Sized + MOCProperties + Iterator<Item=MocRange<T, Q>>
  where
    T: Idx,
    Q: MocQty<T>,
{
  /// Transforms this iterator over ranges into a struct implementing the trait `MOCIterator`.
  fn to_cells_iter(self) -> CellMOCIteratorFromRanges<T, Q, Self> {
    CellMOCIteratorFromRanges::new(self)
  }
}
impl<T, Q, R> RangeMOCIterator<T, Q> for R
  where
    T: Idx,
    Q: MocQty<T>,
    R: Sized + MOCProperties + Iterator<Item=MocRange<T, Q>> {}
*/
pub trait RangeMOCIterator<T: Idx>: Sized + MOCProperties + Iterator<Item=Range<T>> {
  type Qty: MocQty<T>;

  fn cells(self) -> CellMOCIteratorFromRanges<T, Self::Qty, Self> {
    CellMOCIteratorFromRanges::new(self)
  }
}

pub trait RangeMOCIntoIterator<T: Idx>: Sized {
  type Qty: MocQty<T>;
  type IntoRangeMOCIter: RangeMOCIterator<T, Qty=Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter;
}

/// Transforms a `RangeMOCIterator` into a `CellMOCIterator`.
pub struct CellMOCIteratorFromRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  it: R,
  curr: Option<MocRange<T, Q>>,
  depth_max: u8,
  shift_dd: usize,
  range_len_min: T,
  mask: T,
}

impl<T, Q, R> CellMOCIteratorFromRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{

  fn new(mut it: R) -> CellMOCIteratorFromRanges<T, Q, R> {
    let curr = it.next().map(|range| range.into());
    let depth_max = it.depth_max();
    let shift_dd = Q::shift_from_depth_max (depth_max) as usize;
    let range_len_min = T::one() << shift_dd;
    let mask = From::from(Q::LEVEL_MASK << shift_dd);
    CellMOCIteratorFromRanges {
      it,
      curr,
      depth_max,
      shift_dd,
      range_len_min,
      mask,
    }
  }
}

impl<T: Idx, Q: MocQty<T>, R: RangeMOCIterator<T, Qty=Q>> HasMaxDepth for CellMOCIteratorFromRanges<T, Q, R> {
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T: Idx, Q: MocQty<T>, R: RangeMOCIterator<T, Qty=Q>> Sorted for CellMOCIteratorFromRanges<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: RangeMOCIterator<T, Qty=Q>> NonOverlapping for CellMOCIteratorFromRanges<T, Q, R> { }

impl<T, Q, R> Iterator for CellMOCIteratorFromRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{

  type Item = Cell<T>;

  fn next(&mut self) -> Option<Self::Item> {
    if let Some(c) = &mut self.curr {
      let res = c.next_with_knowledge(self.depth_max, self.shift_dd, self.range_len_min, self.mask);
      if res.is_none() {
        self.curr = self.it.next().map(|range| range.into());
        self.next()
      } else {
        res.map(|mc| mc.into())
      }
    } else {
      None
    }
  }
}


// RANGE MOC IMPLEMENTATION

pub struct RangeMOC<T: Idx, Q: MocQty<T>> {
  pub depth_max: u8,
  pub ranges: MocRanges<T, Q>
}
impl<T: Idx, Q: MocQty<T>> RangeMOC<T, Q> {

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
impl<T: Idx, Q: MocQty<T>> Sorted for RangeMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for RangeMOC<T, Q> { }


pub struct RangeMocIter<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: IntoIter<Range<T>>,
  _qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for RangeMocIter<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> Sorted for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> Iterator for RangeMocIter<T, Q> {
  type Item = Range<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for RangeMocIter<T, Q> {
  type Qty = Q;
}
impl<T: Idx, Q: MocQty<T>> RangeMOCIntoIterator<T> for RangeMOC<T, Q> {
  type Qty = Q;
  type IntoRangeMOCIter = RangeMocIter<T, Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter {
    RangeMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.into_iter(),
      _qty: PhantomData
    }
  }
}

pub struct RangeRefMocIter<'a, T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: slice::Iter<'a, Range<T>>,
  _qty: PhantomData<Q>,
}
impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for RangeRefMocIter<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<'a, T: Idx, Q: MocQty<T>> Sorted for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> Iterator for RangeRefMocIter<'a, T, Q> {
  type Item = Range<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(|e| e.clone())
  }
}
impl<'a, T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for RangeRefMocIter<'a, T, Q> {
  type Qty = Q;
}
impl<'a, T: Idx, Q: MocQty<T>> RangeMOCIntoIterator<T> for &'a RangeMOC<T, Q> {
  type Qty = Q;
  type IntoRangeMOCIter = RangeRefMocIter<'a, T, Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter {
    RangeRefMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.iter(),
      _qty: PhantomData
    }
  }
}

/*
/// A very basic, unchecked, Range MOC iterator.
/// We call it lazy because if the stored iterator is the result of a operation,
/// the operation will be performed only as we iterate over the LazyIterator: this
/// is an important point to be aware of when performing benches.
/// For an structure actually performing the operation, see `OwnedRangeMOC`.
pub struct LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  depth_max: u8,
  it: R,
  _t_t: PhantomData<T>,
  _t_q: PhantomData<Q>
}

impl<T, Q, R> LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  pub fn new(depth_max: u8, it: R) -> LazyRangeMOCIter<T, Q, R> {
    LazyRangeMOCIter { depth_max, it, _t_t: PhantomData, _t_q: PhantomData }
  }
}

impl<T, Q, R> HasMaxDepth for LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<T, Q, R> Iterator for LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next()
  }
}
*/




/*
pub struct LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {
  depth_max: u8,
  iter: Iter<'a, HpxRange<H>>
}

impl<'a, H> LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {
  pub fn new(depth_max: u8, iter: Iter<'a, HpxRange<H>>) -> LazyRangeMOCVecIter<'a, H> {
    LazyRangeMOCVecIter{ depth_max, iter}
  }
}

impl<'a, H> HasMaxDepth for LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<'a, H> Iterator for LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {

  type Item = HpxRange<H>;

  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(| HpxRange{ from, to} | HpxRange{ from: *from, to: *to})
  }
}
*/

#[cfg(test)]
mod tests {
  use std::marker::PhantomData;

  use crate::moc::{RangeMOC, RangeMocIter, RangeMOCIntoIterator, CellMOCIteratorFromRanges, RangeMOCIterator};
  use crate::mocranges::MocRanges;
  use crate::ranges::Ranges;
  use crate::mocqty::Hpx;
  use crate::mocell::Cell;
  use crate::hpxranges::HpxUniq2DepthIdxIter;

  #[test]
  fn test_range2cells_1() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![0..5])
    };
    let rit = (&rm).into_range_moc_iter();
    let v1: Vec<Cell<u64>> = rit.cells().collect();
    assert_eq!(v1, vec![
      Cell::new(28, 0),
      Cell::new(29, 4)
    ]);

    let v2: Vec<(i8, u64)> = HpxUniq2DepthIdxIter::new(rm.ranges).collect();
    assert_eq!(v2, vec![
      (28, 0),
      (29, 4)
    ]);
  }

  #[test]
  fn test_range2cells_2() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![2..8])
    };
    let rit = (&rm).into_range_moc_iter();
    let v1: Vec<Cell<u64>> = rit.cells().collect();
    assert_eq!(v1, vec![
      Cell::new(29, 2),
      Cell::new(29, 3),
      Cell::new(28, 1)
    ]);

    let v2: Vec<(i8, u64)> = HpxUniq2DepthIdxIter::new(rm.ranges).collect();
    assert_eq!(v2, vec![
      (28, 1),
      (29, 2),
      (29, 3)
    ]);
  }

  #[test]
  fn test_range2cells_3() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![
        0..5,
        6..59,
        78..6953,
        12458..55587,
        55787..65587
      ])
    };
    let rit = (&rm).into_range_moc_iter();
    let v1: Vec<Cell<u64>> = rit.cells().collect();
    println!("{:?}", v1);

    let v2: Vec<(i8, u64)> = HpxUniq2DepthIdxIter::new(rm.ranges).collect();
    println!("{:?}", v2);

  }

}