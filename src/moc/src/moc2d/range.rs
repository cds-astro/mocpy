
use std::slice;
use std::vec::{IntoIter};

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{
  ZSorted, NonOverlapping,
  RangeMOCIterator, RangeMOCIntoIterator,
  range::{RangeMOC, RangeMocIter, RangeRefMocIter},
  adapters::CellMOCIteratorFromRanges,
};
use crate::moc2d::{
  HasTwoMaxDepth, MOC2Properties, 
  RangeMOC2ElemIt, RangeMOC2Iterator, RangeMOC2IntoIterator,
  CellMOC2ElemIt, CellMOC2Iterator, CellMOC2IntoIterator
};

/// One element of a MOC2 made of Range elements
/// # Info
/// This implementation is not the most memory efficient since it is based on a couple of MOCs,
/// and the `max_depth` of each MOC is thus stored in each element).
/// As an alternative, we could store `MocRanges` insead of `RangeMOC`.
/// We could also replace `Vec<T>` by `Box<[T]>` in the base type `Ranges<T>`
/// TODO later if we run into memory consumption troubles.
#[derive(Debug, Clone)]
pub struct RangeMOC2Elem<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  moc_l: RangeMOC<T, Q>, // or (do not contains the depth): MocRanges<T, Q>,
  moc_r: RangeMOC<U, R>, // or (do not contains the depth): MocRanges<U, R>,
}
impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> RangeMOC2Elem<T, Q, U, R> {
  pub fn new(moc_l: RangeMOC<T, Q>, moc_r: RangeMOC<U, R>) -> Self {
    Self { moc_l , moc_r }
  }
  /// Returns the number of ranges in the first dimension
  pub fn n_ranges_1(&self) -> usize {
    self.moc_l.len()
  }
  /// Returns the number of ranges in the second dimension
  pub fn n_ranges_2(&self) -> usize {
    self.moc_r.len()
  }
  /// Returns the number of ranges in both quantities
  pub fn n_ranges(&self) -> u64 {
    self.n_ranges_1() as u64 + self.n_ranges_2() as u64
  }
  
  pub fn mocs(self) -> (RangeMOC<T, Q>, RangeMOC<U, R>) {
    (self.moc_l, self.moc_r)
  }
}
impl<T, Q, U, R> RangeMOC2ElemIt<T, Q, U, R> for RangeMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = RangeMocIter<T, Q>;
  type It2 = RangeMocIter<U, R>;
  fn range_mocs_it(self) -> (Self::It1,  Self::It2)  {
    (self.moc_l.into_range_moc_iter(), self.moc_r.into_range_moc_iter())
  }
}
impl<'a, T, Q, U, R> RangeMOC2ElemIt<T, Q, U, R> for &'a RangeMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = RangeRefMocIter<'a, T, Q>;
  type It2 = RangeRefMocIter<'a, U, R>;
  fn range_mocs_it(self) -> (Self::It1,  Self::It2)  {
    ((&self.moc_l).into_range_moc_iter(), (&self.moc_r).into_range_moc_iter())
  }
}


impl<T, Q, U, R> CellMOC2ElemIt<T, Q, U, R> for RangeMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = CellMOCIteratorFromRanges<T, Q, RangeMocIter<T, Q>> ;
  type It2 = CellMOCIteratorFromRanges<U, R, RangeMocIter<U, R>>;
  fn cell_mocs_it(self) -> (Self::It1,  Self::It2)  {
    (self.moc_l.into_range_moc_iter().cells(), self.moc_r.into_range_moc_iter().cells())
  }
}
impl<'a, T, Q, U, R> CellMOC2ElemIt<T, Q, U, R> for &'a RangeMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = CellMOCIteratorFromRanges<T, Q, RangeRefMocIter<'a, T, Q>> ;
  type It2 = CellMOCIteratorFromRanges<U, R, RangeRefMocIter<'a, U, R>>;
  fn cell_mocs_it(self) -> (Self::It1,  Self::It2)  {
    ((&self.moc_l).into_range_moc_iter().cells(), (&self.moc_r).into_range_moc_iter().cells())
  }
}


/// A MOC2 made of Range elements
pub struct RangeMOC2<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  depth_max_l: u8,
  depth_max_r: u8, // not in vmoc. Really useful?
  elems: Vec<RangeMOC2Elem<T, Q, U, R>>,
}
impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>  RangeMOC2<T, Q, U, R> {
  pub fn new(depth_max_l: u8, depth_max_r: u8, elems: Vec<RangeMOC2Elem<T, Q, U, R>>) -> Self {
    Self { depth_max_l, depth_max_r, elems }
  }
  /// The total number of ranges in both dimensions
  pub fn compute_n_ranges(&self) -> u64 {
    self.elems.iter().map(|e| e.n_ranges()).sum()
  }
}
impl<T, Q, U, R> HasTwoMaxDepth for RangeMOC2<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  fn depth_max_1(&self) -> u8 {
    self.depth_max_l
  }
  fn depth_max_2(&self) -> u8 {
    self.depth_max_r
  }
}
impl<T, Q, U, R> ZSorted for RangeMOC2<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> NonOverlapping for RangeMOC2<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> MOC2Properties for RangeMOC2<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}

/// Iterator taking the ownership of a MOC2 made of Range elements
pub struct RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  depth_max_l: u8,
  depth_max_r: u8,
  iter: IntoIter<RangeMOC2Elem<T, Q, U, R>>
}
impl<T, Q, U, R> HasTwoMaxDepth for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  fn depth_max_1(&self) -> u8 {
    self.depth_max_l
  }
  fn depth_max_2(&self) -> u8 {
    self.depth_max_r
  }
}
impl<T, Q, U, R> ZSorted for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> NonOverlapping for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> MOC2Properties for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> Iterator for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type Item = RangeMOC2Elem<T, Q, U, R>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}


impl<T, Q, U, R> RangeMOC2Iterator<
  T, Q, RangeMocIter<T, Q>,
  U, R, RangeMocIter<U, R>,
  RangeMOC2Elem<T, Q, U, R>
> for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> { }

impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>
RangeMOC2IntoIterator<
  T, Q, RangeMocIter<T, Q>,
  U, R, RangeMocIter<U, R>,
  RangeMOC2Elem<T, Q, U, R>
> for RangeMOC2<T, Q, U, R> {
  
  type IntoRangeMOC2Iter = RangeMoc2Iter<T, Q, U, R>;
  
  fn into_range_moc2_iter(self) -> Self::IntoRangeMOC2Iter {
    RangeMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.into_iter()
    }
  }
}


impl<T, Q, U, R> CellMOC2Iterator<
  T, Q, CellMOCIteratorFromRanges<T, Q, RangeMocIter<T, Q>>,
  U, R, CellMOCIteratorFromRanges<U, R, RangeMocIter<U, R>>,
  RangeMOC2Elem<T, Q, U, R>
> for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> { }

impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>
CellMOC2IntoIterator<
  T, Q, CellMOCIteratorFromRanges<T, Q, RangeMocIter<T, Q>>,
  U, R, CellMOCIteratorFromRanges<U, R, RangeMocIter<U, R>>,
  RangeMOC2Elem<T, Q, U, R>
> for RangeMOC2<T, Q, U, R> {

  type IntoCellMOC2Iter = RangeMoc2Iter<T, Q, U, R>;

  fn into_cell_moc2_iter(self) -> Self::IntoCellMOC2Iter {
    RangeMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.into_iter()
    }
  }
}


/*impl<T, Q, U, R> CellMOC2Iterator<
  T, Q, CellMOCIteratorFromRanges<T, Q, RangeMocIter<T, Q>>,
  U, R, CellMOCIteratorFromRanges<U, R, RangeMocIter<U, R>>,
  RangeMOC2Elem<T, Q, U, R>
> for RangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}

impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> 
CellMOC2IntoIterator<
  T, Q, CellMOCIteratorFromRanges<T, Q, RangeMocIter<T, Q>>,
  U, U, CellMOCIteratorFromRanges<U, R, RangeMocIter<U, R>>,
  RangeMOC2ElemIt<
    T, Q, U, R, 
    It1=CellMOCIteratorFromRanges<T, Q, RangeMocIter<T, Q>>,
    It2=CellMOCIteratorFromRanges<U, R, RangeMocIter<U, R>>
  >
> for RangeMOC2<T, Q, U, R> {
  type IntoCellMOC2Iter = RangeMoc2Iter<T, Q, U, R>;

  fn into_cell_moc2_iter(self) -> Self::IntoCellMOC2Iter {
    RangeMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.into_iter()
    }
  }
}*/

/// Iterator borrowing a MOC2 made of Range elements
pub struct RangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  depth_max_l: u8,
  depth_max_r: u8,
  iter: slice::Iter<'a, RangeMOC2Elem<T, Q, U, R>>
}
impl<'a, T, Q, U, R> HasTwoMaxDepth for RangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  fn depth_max_1(&self) -> u8 {
    self.depth_max_l
  }
  fn depth_max_2(&self) -> u8 {
    self.depth_max_r
  }
}
impl<'a, T, Q, U, R> ZSorted for RangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> NonOverlapping for RangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> MOC2Properties for RangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> Iterator for RangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type Item = &'a RangeMOC2Elem<T, Q, U, R>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}

impl<'a, T, Q, U, R> RangeMOC2Iterator<
  T, Q, RangeRefMocIter<'a, T, Q>,
  U, R, RangeRefMocIter<'a, U, R>,
  &'a RangeMOC2Elem<T, Q, U, R>
> for RangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}

impl<'a, T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>
RangeMOC2IntoIterator<
  T, Q, RangeRefMocIter<'a, T, Q>,
  U, R, RangeRefMocIter<'a, U, R>,
  &'a RangeMOC2Elem<T, Q, U, R>
> for &'a RangeMOC2<T, Q, U, R> {
  type IntoRangeMOC2Iter = RangeRefMoc2Iter<'a, T, Q, U, R>;
  fn into_range_moc2_iter(self) -> Self::IntoRangeMOC2Iter {
    RangeRefMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.iter()
    }
  }
}
