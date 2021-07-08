
use std::slice;
use std::vec::IntoIter;
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::elem::cell::Cell;
use crate::elemset::cell::MocCells;
use crate::moc::{
  HasMaxDepth, ZSorted, NonOverlapping, MOCProperties,
  CellMOCIterator, CellMOCIntoIterator
};

/// A MOC made of (ordered and non-overlaping) cells.
/// This is used as the result of a MOC JSON deserialization of a MOC.
#[derive(Debug)]
pub struct CellMOC<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  cells: MocCells<T, Q>
}
impl<T: Idx, Q: MocQty<T>> CellMOC<T, Q> {
  pub fn new(depth_max: u8, cells: MocCells<T, Q>) -> Self {
    Self { depth_max, cells }
  }
  pub fn len(&self) -> usize {
    self.cells.cells().len()
  }
  pub fn is_empty(&self) -> bool { self.len() == 0 }

}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellMOC<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for CellMOC<T, Q> { }

/// Iterator taking the ownership of the `CellMOC` it iterates over.
pub struct CellMocIter<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  last: Option<Cell<T>>,
  iter: IntoIter<Cell<T>>,
  _qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellMocIter<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for CellMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> Iterator for CellMocIter<T, Q> {
  type Item = Cell<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<T: Idx, Q: MocQty<T>> CellMOCIterator<T> for CellMocIter<T, Q> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Cell<T>> {
    self.last.as_ref()
  }
}
impl<T: Idx, Q: MocQty<T>> CellMOCIntoIterator<T> for CellMOC<T, Q> {
  type Qty = Q;
  type IntoCellMOCIter = CellMocIter<T, Self::Qty>;

  fn into_cell_moc_iter(self) -> Self::IntoCellMOCIter {
    let l = self.cells.0.0.len();
    let last: Option<Cell<T>> = if l > 0 {
      Some(self.cells.0.0[l - 1])
    } else {
      None
    };
    CellMocIter {
      depth_max: self.depth_max,
      last,
      iter: self.cells.0.0.into_vec().into_iter(),
      _qty: PhantomData
    }
  }
}

/// Iterator borrowing the `CellMOC` it iterates over.
pub struct CellRefMocIter<'a, T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  last: Option<Cell<T>>,
  iter: slice::Iter<'a, Cell<T>>,
  _qty: PhantomData<Q>,
}
impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for CellRefMocIter<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<'a, T: Idx, Q: MocQty<T>> ZSorted for CellRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for CellRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for CellRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> Iterator for CellRefMocIter<'a, T, Q> {
  type Item = Cell<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().cloned()
  }
}
impl<'a, T: Idx, Q: MocQty<T>> CellMOCIterator<T> for CellRefMocIter<'a, T, Q> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Cell<T>> {
    self.last.as_ref()
  }
}
impl<'a, T: Idx, Q: MocQty<T>> CellMOCIntoIterator<T> for &'a CellMOC<T, Q> {
  type Qty = Q;
  type IntoCellMOCIter = CellRefMocIter<'a, T, Self::Qty>;

  fn into_cell_moc_iter(self) -> Self::IntoCellMOCIter {
    let l = self.cells.0.0.len();
    let last: Option<Cell<T>> = if l > 0 {
      Some(self.cells.0.0[l - 1])
    } else {
      None
    };
    CellRefMocIter {
      depth_max: self.depth_max,
      last,
      iter: self.cells.0.0.iter(),
      _qty: PhantomData
    }
  }
}
