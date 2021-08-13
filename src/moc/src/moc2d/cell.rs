
use std::slice;
use std::vec::{IntoIter};

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{
  ZSorted, NonOverlapping, CellMOCIntoIterator,
  cell::{CellMOC, CellMocIter, CellRefMocIter}
};
use crate::moc2d::{
  HasTwoMaxDepth, MOC2Properties,
  CellMOC2ElemIt, CellMOC2Iterator, CellMOC2IntoIterator
};

/// One element of a MOC2 made of Cell elements.
/// # Info
/// It would be more memory efficient to store internally `MocCells` objects
/// instead of `CellMOC` objects (because of the `max_depth`).
/// But then we need a way to provide the `depth_max` when serializing the elements...
/// For memory efficiency, see `RangeMOC2`.
pub struct CellMOC2Elem<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  moc_l: CellMOC<T, Q>, // or (do not contains the depth): MocCells<T, Q>,
  moc_r: CellMOC<U, R>, // or (do not contains the depth): MocCells<U, R>,
}
impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>  CellMOC2Elem<T, Q, U, R> {
  pub fn new(moc_l: CellMOC<T, Q>, moc_r: CellMOC<U, R>) -> Self {
    Self { moc_l , moc_r }
  }

  pub fn mocs(self) -> (CellMOC<T, Q>, CellMOC<U, R>) {
    (self.moc_l, self.moc_r)
  }
}
impl<T, Q, U, R> CellMOC2ElemIt<T, Q, U, R> for CellMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = CellMocIter<T, Q>;
  type It2 = CellMocIter<U, R>;
  fn cell_mocs_it(self) -> (Self::It1,  Self::It2)  {
    (self.moc_l.into_cell_moc_iter(), self.moc_r.into_cell_moc_iter())
  }
}
impl<'a, T, Q, U, R> CellMOC2ElemIt<T, Q, U, R> for &'a CellMOC2Elem<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type It1 = CellRefMocIter<'a, T, Q>;
  type It2 = CellRefMocIter<'a, U, R>;
  fn cell_mocs_it(self) -> (Self::It1,  Self::It2)  {
    ((&self.moc_l).into_cell_moc_iter(), (&self.moc_r).into_cell_moc_iter())
  }
}

/// A MOC2 made of Cell elements
pub struct CellMOC2<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  depth_max_l: u8,
  depth_max_r: u8, // not in vmoc. Really useful?
  elems: Vec<CellMOC2Elem<T, Q, U, R>>,
}
impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>  CellMOC2<T, Q, U, R> {
  pub fn new(depth_max_l: u8, depth_max_r: u8, elems: Vec<CellMOC2Elem<T, Q, U, R>>) -> Self {
    Self { depth_max_l, depth_max_r, elems }
  }
  pub fn n_entries(&self) -> usize {
    self.elems.len()
  }
}

impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> HasTwoMaxDepth for CellMOC2<T, Q, U, R> {
  fn depth_max_1(&self) -> u8 {
    self.depth_max_l
  }
  fn depth_max_2(&self) -> u8 {
    self.depth_max_r
  }
}

/// Iterator taking the ownership of a MOC2 made of Cell elements
pub struct CellMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  depth_max_l: u8,
  depth_max_r: u8,
  iter: IntoIter<CellMOC2Elem<T, Q, U, R>>
}
impl<T, Q, U, R> HasTwoMaxDepth for CellMoc2Iter<T, Q, U, R>
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
impl<T, Q, U, R> ZSorted for CellMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> NonOverlapping for CellMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> MOC2Properties for CellMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<T, Q, U, R> Iterator for CellMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type Item = CellMOC2Elem<T, Q, U, R>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<T, Q, U, R> CellMOC2Iterator<
  T, Q, CellMocIter<T, Q>,
  U, R, CellMocIter<U, R>,
  CellMOC2Elem<T, Q, U, R>
> for CellMoc2Iter<T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}

impl<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>
CellMOC2IntoIterator<
  T, Q, CellMocIter<T, Q>,
  U, R, CellMocIter<U, R>,
  CellMOC2Elem<T, Q, U, R>
> for CellMOC2<T, Q, U, R> {
  type IntoCellMOC2Iter = CellMoc2Iter<T, Q, U, R>;

  fn into_cell_moc2_iter(self) -> Self::IntoCellMOC2Iter {
    CellMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.into_iter()
    }
  }
}

/// Iterator borrowing a MOC2 made of Cell elements
pub struct CellRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  depth_max_l: u8,
  depth_max_r: u8,
  iter: slice::Iter<'a, CellMOC2Elem<T, Q, U, R>>
}
impl<'a, T, Q, U, R> HasTwoMaxDepth for CellRefMoc2Iter<'a, T, Q, U, R>
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
impl<'a, T, Q, U, R> ZSorted for CellRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> NonOverlapping for CellRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> MOC2Properties for CellRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T, Q, U, R> Iterator for CellRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>
{
  type Item = &'a CellMOC2Elem<T, Q, U, R>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<'a, T, Q, U, R> CellMOC2Iterator<
  T, Q, CellRefMocIter<'a, T, Q>,
  U, R, CellRefMocIter<'a, U, R>,
  &'a CellMOC2Elem<T, Q, U, R>
> for CellRefMoc2Iter<'a, T, Q, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U> {}
impl<'a, T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>>
CellMOC2IntoIterator<
  T, Q, CellRefMocIter<'a, T, Q>,
  U, R, CellRefMocIter<'a, U, R>,
  &'a CellMOC2Elem<T, Q, U, R>
> for &'a CellMOC2<T, Q, U, R> {
  type IntoCellMOC2Iter = CellRefMoc2Iter<'a, T, Q, U, R>;

  fn into_cell_moc2_iter(self) -> Self::IntoCellMOC2Iter {
    CellRefMoc2Iter {
      depth_max_l: self.depth_max_l,
      depth_max_r: self.depth_max_r,
      iter: self.elems.iter()
    }
  }
}
