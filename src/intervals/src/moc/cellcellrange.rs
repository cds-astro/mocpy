
use std::slice;
use std::vec::{IntoIter};
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, CellOrCellRangeMOCIterator, CellOrCellRangeMOCIntoIterator};
use crate::mocells::{MocCellOrCellRanges, CellOrCellRanges};
use crate::mocell::CellOrCellRange;

/// A MOC made of a mix of (ordered and non-overlaping) cells and cells range
/// (a cell range is a range at the cell depth, while regular ranges are at the
/// largest possible depth and depends on type  `T`).
/// This is used as the result of a MOC ASCII deserialization.
pub struct CellOrCellRangeMOC<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  ranges: MocCellOrCellRanges<T, Q>
}
impl<T: Idx, Q: MocQty<T>> CellOrCellRangeMOC<T, Q> {
  pub fn new(depth_max: u8, ranges: MocCellOrCellRanges<T, Q>) -> Self {
    Self { depth_max, ranges }
  }
  pub fn into_cellcellrange_moc_iter(self) -> CellOrCellRangeMocIter<T, Q> {
    CellOrCellRangeMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.into_iter(),
      _qty: PhantomData
    }
  }
  pub fn elems(self) -> CellOrCellRanges<T> {
    self.ranges.0
  }
  pub fn moc_elems(self) -> MocCellOrCellRanges<T, Q> {
    self.ranges
  }
}

impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellOrCellRangeMOC<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellOrCellRangeMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellOrCellRangeMOC<T, Q> { }

/// Iterator taking the ownership of the `CellOrCellRangeMOC` it iterates over.
pub struct CellOrCellRangeMocIter<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: IntoIter<CellOrCellRange<T>>,
  _qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellOrCellRangeMocIter<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellOrCellRangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellOrCellRangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for CellOrCellRangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> Iterator for CellOrCellRangeMocIter<T, Q> {
  type Item = CellOrCellRange<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIterator<T> for CellOrCellRangeMocIter<T, Q> {
  type Qty = Q;
}
impl<T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIntoIterator<T> for CellOrCellRangeMOC<T, Q> {
  type Qty = Q;
  type IntoCellOrCellRangeMOCIter = CellOrCellRangeMocIter<T, Self::Qty>;

  fn into_cellcellrange_moc_iter(self) -> Self::IntoCellOrCellRangeMOCIter {
    CellOrCellRangeMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.into_iter(),
      _qty: PhantomData
    }
  }
}

/// Iterator borrowing the `CellOrCellRangeMOC` it iterates over.
pub struct CellOrCellRangeRefMocIter<'a, T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: slice::Iter<'a, CellOrCellRange<T>>,
  _qty: PhantomData<Q>,
}
impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for CellOrCellRangeRefMocIter<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<'a, T: Idx, Q: MocQty<T>> ZSorted for CellOrCellRangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for CellOrCellRangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for CellOrCellRangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> Iterator for CellOrCellRangeRefMocIter<'a, T, Q> {
  type Item = CellOrCellRange<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(|e| e.clone())
  }
  // Declaring size_hint, a 'collect' can directly allocate the right number of elements
  fn size_hint(&self) -> (usize, Option<usize>) {
    self.iter.size_hint()
  }
}
impl<'a, T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIterator<T> for CellOrCellRangeRefMocIter<'a, T, Q> {
  type Qty = Q;
}
impl<'a, T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIntoIterator<T> for &'a CellOrCellRangeMOC<T, Q> {
  type Qty = Q;
  type IntoCellOrCellRangeMOCIter = CellOrCellRangeRefMocIter<'a, T, Self::Qty>;

  fn into_cellcellrange_moc_iter(self) -> Self::IntoCellOrCellRangeMOCIter {
    CellOrCellRangeRefMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.iter(),
      _qty: PhantomData
    }
  }
}
