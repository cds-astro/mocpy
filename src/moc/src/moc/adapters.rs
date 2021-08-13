
use std::ops::Range;

use crate::idx::Idx;
use super::{
  HasMaxDepth,
  RangeMOCIterator
};
use crate::qty::MocQty;
use crate::elem::{
  cell::Cell,
  cellrange::CellRange,
  cellcellrange::CellOrCellRange,
  range::MocRange,
};
use crate::moc::{CellMOCIterator, ZSorted, NonOverlapping, MOCProperties, CellOrCellRangeMOCIterator};

/// Transforms a `RangeMOCIterator` into a `CellMOCIterator`.
pub struct CellMOCIteratorFromRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  it: R,
  last: Option<Cell<T>>,
  curr: Option<MocRange<T, Q>>,
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

  pub fn new(mut it: R) -> CellMOCIteratorFromRanges<T, Q, R> {
    let last = it.peek_last().and_then(|r| {
      let moc_range: MocRange<T, Q> = r.into();
      // TODO: Make a method retrieving directly the last elem instead of iterating
      moc_range.last().map(|mc| mc.into())
    });
    let curr = it.next().map(|range| range.into());
    let shift_dd = Q::shift_from_depth_max (it.depth_max()) as usize;
    let range_len_min = T::one() << shift_dd;
    let mut mask: T = From::from(Q::LEVEL_MASK);
    mask = mask.unsigned_shl(shift_dd as u32);
    CellMOCIteratorFromRanges {
      it,
      last,
      curr,
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
impl<T: Idx, Q: MocQty<T>, R: RangeMOCIterator<T, Qty=Q>> ZSorted for CellMOCIteratorFromRanges<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: RangeMOCIterator<T, Qty=Q>> NonOverlapping for CellMOCIteratorFromRanges<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: RangeMOCIterator<T, Qty=Q>> MOCProperties for CellMOCIteratorFromRanges<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: RangeMOCIterator<T, Qty=Q>> CellMOCIterator<T> for CellMOCIteratorFromRanges<T, Q, R> {
  type Qty = Q;
  
  fn peek_last(&self) -> Option<&Cell<T>> {
    self.last.as_ref()
  } 
}
impl<T, Q, R> Iterator for CellMOCIteratorFromRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{

  type Item = Cell<T>;

  fn next(&mut self) -> Option<Self::Item> {
    if let Some(c) = &mut self.curr {
      let res = c.next_cell_with_knowledge(self.it.depth_max(), self.shift_dd, self.range_len_min, self.mask);
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

/// Transforms a `CellMOCIterator` into a `RangeMOCIterator`.
pub struct RangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  last: Option<Range<T>>,
  it: R,
  curr: Option<Range<T>>,
}
impl<T, Q, R> RangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  pub fn new(mut it: R, last: Option<Range<T>>) -> RangeMOCIteratorFromCells<T, Q, R> {
    let curr: Option<Range<T>> = it.next().map(|e| MocRange::<T, Q>::from(e).0);
    RangeMOCIteratorFromCells {
      last,
      it,
      curr,
    }
  }
}

impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> HasMaxDepth for RangeMOCIteratorFromCells<T, Q, R> {
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> ZSorted for RangeMOCIteratorFromCells<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> NonOverlapping for RangeMOCIteratorFromCells<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> MOCProperties for RangeMOCIteratorFromCells<T, Q, R> { }

impl<T, Q, R> Iterator for RangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    match &mut self.curr {
      Some(Range { start: _, end: lend }) => {
        let mut next = self.it.next().map(|e| MocRange::<T, Q>::from(e).0);
        while let Some(Range { start: rstart, end: rend }) = next {
          if rstart <= *lend {
            *lend = rend;
            next = self.it.next().map(|e| MocRange::<T, Q>::from(e).0);
          } else {
            break;
          }
        }
        std::mem::replace(&mut self.curr, next)
      },
      None => None,
    }
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let (low, upp) = self.it.size_hint();
    if low > 0 {
      (1, upp)
    } else {
      (0, upp)
    }
  }
}
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> RangeMOCIterator<T> for RangeMOCIteratorFromCells<T, Q, R> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    self.last.as_ref()
  }
}


/// Transforms a `CellOrCellRangeMOCIterator` into a `RangeMOCIterator`.
pub struct RangeMOCIteratorFromCellOrCellRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellOrCellRangeMOCIterator<T, Qty=Q>
{
  last: Option<Range<T>>,
  it: R,
  curr: Option<Range<T>>,
}
impl<T, Q, R> RangeMOCIteratorFromCellOrCellRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellOrCellRangeMOCIterator<T, Qty=Q>
{
  pub fn new(mut it: R, last: Option<Range<T>>) -> RangeMOCIteratorFromCellOrCellRanges<T, Q, R> {
    let curr: Option<Range<T>> = it.next().map(|e| MocRange::<T, Q>::from(e).0);
    RangeMOCIteratorFromCellOrCellRanges {
      last,
      it,
      curr,
    }
  }
}

impl<T: Idx, Q: MocQty<T>, R: CellOrCellRangeMOCIterator<T, Qty=Q>> HasMaxDepth for RangeMOCIteratorFromCellOrCellRanges<T, Q, R> {
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T: Idx, Q: MocQty<T>, R: CellOrCellRangeMOCIterator<T, Qty=Q>> ZSorted for RangeMOCIteratorFromCellOrCellRanges<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: CellOrCellRangeMOCIterator<T, Qty=Q>> NonOverlapping for RangeMOCIteratorFromCellOrCellRanges<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: CellOrCellRangeMOCIterator<T, Qty=Q>> MOCProperties for RangeMOCIteratorFromCellOrCellRanges<T, Q, R> { }

impl<T, Q, R> Iterator for RangeMOCIteratorFromCellOrCellRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellOrCellRangeMOCIterator<T, Qty=Q>
{
  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    match &mut self.curr {
      Some(Range { start: _, end: lend }) => {
        let mut next = self.it.next().map(|e| MocRange::<T, Q>::from(e).0);
        while let Some(Range { start: rstart, end: rend }) = next {
          if rstart <= *lend {
            *lend = rend;
            next = self.it.next().map(|e| MocRange::<T, Q>::from(e).0);
          } else {
            break;
          }
        }
        std::mem::replace(&mut self.curr, next)
      },
      None => None,
    }
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let (low, upp) = self.it.size_hint();
    if low > 0 {
      (1, upp)
    } else {
      (0, upp)
    }
  }
}
impl<T: Idx, Q: MocQty<T>, R: CellOrCellRangeMOCIterator<T, Qty=Q>> RangeMOCIterator<T> for RangeMOCIteratorFromCellOrCellRanges<T, Q, R> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    self.last.as_ref()
  }
}


/// Transforms a `CellMOCIterator` into a `CellOrCellRangeMOCIterator`.
pub struct CellOrCellRangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  it: R,
  last: Option<CellOrCellRange<T>>,
  curr: Option<Cell<T>>,
}

impl<T, Q, R> CellOrCellRangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  pub fn new(mut it: R) -> CellOrCellRangeMOCIteratorFromCells<T, Q, R> {
    let last = it.peek_last().and_then(|r| {
      let moc_range: MocRange<T, Q> = r.into();
      // TODO: Make a method retrieving directly the last elem instead of iterating
      moc_range.last().map(|mc| CellOrCellRange::Cell(mc.into()))
    });
    let curr = it.next();
    CellOrCellRangeMOCIteratorFromCells {
      it,
      last,
      curr,
    }
  }
}

impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> HasMaxDepth for CellOrCellRangeMOCIteratorFromCells<T, Q, R> {
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> ZSorted for CellOrCellRangeMOCIteratorFromCells<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> NonOverlapping for CellOrCellRangeMOCIteratorFromCells<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> MOCProperties for CellOrCellRangeMOCIteratorFromCells<T, Q, R> { }

impl<T, Q, R> Iterator for CellOrCellRangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  type Item = CellOrCellRange<T>;

  fn next(&mut self) -> Option<Self::Item> {
    if let Some(Cell{ depth: left_d, idx: left_i}) = &self.curr {
      let mut n = T::one();
      let mut next = self.it.next();
      while let Some(Cell { depth: right_d, idx: right_i}) = &next {
        if *left_d == *right_d && *left_i + n == *right_i {
          n += T::one();
          next = self.it.next();
        } else {
          break;
        }
      }
      std::mem::replace(&mut self.curr, next)
        .map(move |c| if n == T::one() {
          CellOrCellRange::Cell(c)
        } else {
          CellOrCellRange::CellRange(CellRange::new(c.depth, c.idx, c.idx + n))
        })
    } else {
      None
    }
  }
  // We do not declare a size_hint so far because we use it only in streaming mode when producing
  // ASCII output.
}
impl<T: Idx, Q: MocQty<T>, R: CellMOCIterator<T, Qty=Q>> CellOrCellRangeMOCIterator<T> for CellOrCellRangeMOCIteratorFromCells<T, Q, R> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&CellOrCellRange<T>> {
    self.last.as_ref()
  }
}
