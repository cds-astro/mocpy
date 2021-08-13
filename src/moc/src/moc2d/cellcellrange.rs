
use std::slice;
use std::vec::{IntoIter};

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{
  ZSorted, NonOverlapping, CellOrCellRangeMOCIntoIterator,
  cellcellrange::{CellOrCellRangeMOC, CellOrCellRangeMocIter, CellOrCellRangeRefMocIter}
};
use crate::moc2d::{
  HasTwoMaxDepth, MOC2Properties,
  CellOrCellRangeMOC2ElemIt, CellOrCellRangeMOC2Iterator, CellOrCellRangeMOC2IntoIterator
};

/// One element of a MOC2 made of CellOrCellRange elements.
/// # Info
/// It would be more memory efficient to store internally `MocCellOrCellRanges` objects
/// instead of `CellOrCellRangeMOC` objects (because of the `max_depth`).
/// But then we need a way to provide the `depth_max` when serializing the elements...
/// For memory efficiency, see `RangeMOC2`.
pub struct CellOrCellRangeMOC2Elem<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  moc_l: CellOrCellRangeMOC<T, Q>, // or (do not contains the depth): MocCellOrCellRanges<T, Q>,
  moc_r: CellOrCellRangeMOC<U, R>, // or (do not contains the depth): MocCellOrCellRanges<U, R>,
}
impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>  CellOrCellRangeMOC2Elem<T, Q, U, R> {
  pub fn new(moc_l: CellOrCellRangeMOC<T, Q>, moc_r: CellOrCellRangeMOC<U, R>) -> Self {
    Self { moc_l , moc_r }
  }

  pub fn mocs(self) -> (CellOrCellRangeMOC<T, Q>, CellOrCellRangeMOC<U, R>) {
    (self.moc_l, self.moc_r)
  }
}
impl<T, Q, U, R> CellOrCellRangeMOC2ElemIt<T, Q, U, R> for CellOrCellRangeMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = CellOrCellRangeMocIter<T, Q>;
  type It2 = CellOrCellRangeMocIter<U, R>;
  fn cellcellrange_mocs_it(self) -> (Self::It1,  Self::It2)  {
    (self.moc_l.into_cellcellrange_moc_iter(), self.moc_r.into_cellcellrange_moc_iter())
  }
}
impl<'a, T, Q, U, R> CellOrCellRangeMOC2ElemIt<T, Q, U, R> for &'a CellOrCellRangeMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = CellOrCellRangeRefMocIter<'a, T, Q>;
  type It2 = CellOrCellRangeRefMocIter<'a, U, R>;
  fn cellcellrange_mocs_it(self) -> (Self::It1,  Self::It2)  {
    ((&self.moc_l).into_cellcellrange_moc_iter(), (&self.moc_r).into_cellcellrange_moc_iter())
  }
}

/// A MOC2 made of CellOrCellRange elements
pub struct CellOrCellRangeMOC2<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  depth_max_l: u8,
  depth_max_r: u8, // not in vmoc. Really useful?
  elems: Vec<CellOrCellRangeMOC2Elem<T, Q, U, R>>,
}
impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> CellOrCellRangeMOC2<T, Q, U, R> {
  pub fn new(depth_max_l: u8, depth_max_r: u8, elems: Vec<CellOrCellRangeMOC2Elem<T, Q, U, R>>) -> Self {
    Self { depth_max_l, depth_max_r, elems }
  }
}
impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> HasTwoMaxDepth for CellOrCellRangeMOC2<T, Q, U, R> {
  fn depth_max_1(&self) -> u8 {
    self.depth_max_l
  }
  fn depth_max_2(&self) -> u8 {
    self.depth_max_r
  }
}



/// Iterator taking the ownership of a MOC2 made of CellOrCellRange elements
pub struct CellOrCellRangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  depth_max_l: u8,
  depth_max_r: u8,
  iter: IntoIter<CellOrCellRangeMOC2Elem<T, Q, U, R>>
}
impl<T, Q, U, R> HasTwoMaxDepth for CellOrCellRangeMoc2Iter<T, Q, U, R>
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
impl<T, Q, U, R> ZSorted for CellOrCellRangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> NonOverlapping for CellOrCellRangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> MOC2Properties for CellOrCellRangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> Iterator for CellOrCellRangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type Item = CellOrCellRangeMOC2Elem<T, Q, U, R>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<T, Q, U, R> CellOrCellRangeMOC2Iterator<
  T, Q, CellOrCellRangeMocIter<T, Q>,
  U, R, CellOrCellRangeMocIter<U, R>,
  CellOrCellRangeMOC2Elem<T, Q, U, R>
> for CellOrCellRangeMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}

impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>
CellOrCellRangeMOC2IntoIterator<
  T, Q, CellOrCellRangeMocIter<T, Q>,
  U, R, CellOrCellRangeMocIter<U, R>,
  CellOrCellRangeMOC2Elem<T, Q, U, R>
> for CellOrCellRangeMOC2<T, Q, U, R> {
  type IntoCellOrCellRangeMOC2Iter = CellOrCellRangeMoc2Iter<T, Q, U, R>;

  fn into_cellcellrange_moc2_iter(self) -> Self::IntoCellOrCellRangeMOC2Iter {
    CellOrCellRangeMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.into_iter()
    }
  }
}

/// Iterator borrowing a MOC2 made of CellOrCellRange elements
pub struct CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  depth_max_l: u8,
  depth_max_r: u8,
  iter: slice::Iter<'a, CellOrCellRangeMOC2Elem<T, Q, U, R>>
}
impl<'a, T, Q, U, R> HasTwoMaxDepth for CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>
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
impl<'a, T, Q, U, R> ZSorted for CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> NonOverlapping for CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> MOC2Properties for CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> Iterator for CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type Item = &'a CellOrCellRangeMOC2Elem<T, Q, U, R>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<'a, T, Q, U, R> CellOrCellRangeMOC2Iterator<
  T, Q, CellOrCellRangeRefMocIter<'a, T, Q>,
  U, R, CellOrCellRangeRefMocIter<'a, U, R>,
  &'a CellOrCellRangeMOC2Elem<T, Q, U, R>
> for CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>
CellOrCellRangeMOC2IntoIterator<
  T, Q, CellOrCellRangeRefMocIter<'a, T, Q>,
  U, R, CellOrCellRangeRefMocIter<'a, U, R>,
  &'a CellOrCellRangeMOC2Elem<T, Q, U, R>
> for &'a CellOrCellRangeMOC2<T, Q, U, R> {
  type IntoCellOrCellRangeMOC2Iter = CellOrCellRangeRefMoc2Iter<'a, T, Q, U, R>;

  fn into_cellcellrange_moc2_iter(self) -> Self::IntoCellOrCellRangeMOC2Iter {
    CellOrCellRangeRefMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.iter()
    }
  }
}
