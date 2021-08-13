
use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};

/// Performs a logical `NOT` on the input iterator of ranges (i.e. returns its complement).
pub fn not<T, Q, I>(it: I) -> NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>
{
  NotRangeIter::new(it)
}

/// Performs a logical `NOT` on the input iterator of ranges (i.e. returns its complement).
pub struct NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  it: I,
  curr: Option<Range<T>>,
  start: T,
  n_cells_max: T,
}

impl<T, Q, I> NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{

  pub fn new(mut it: I) -> NotRangeIter<T, Q, I>  {
    let n_cells_max = Q::n_cells_max();
    match it.next() {
      Some(range) =>
        if range.start == T::zero() {
          if range.end == n_cells_max {
            // input = all domain => Not is empty
            NotRangeIter {
              it,
              curr: None,
              start: n_cells_max,
              n_cells_max
            }
          } else {
            // look at the second range
            match it.next() {
              Some(range2) =>
                NotRangeIter {
                  it,
                  curr: Some(range.end..range2.start),
                  start: range2.end,
                  n_cells_max
                },
              None =>
                NotRangeIter {
                  it,
                  curr: Some(range.end..n_cells_max),
                  start: n_cells_max,
                  n_cells_max
                }
            }
          }
        } else {
          // Case in which we have to add a cell before the first input cell
          NotRangeIter {
            it,
            curr: Some(T::zero()..range.start),
            start: range.end,
            n_cells_max
          }
        },
      None => {
        // The input is empty => Not is the full domain
        NotRangeIter {
          it,
          curr: Some(T::zero()..n_cells_max),
          start: n_cells_max,
          n_cells_max
        }
      },
    }
  }
}

impl<T, Q, I> HasMaxDepth for  NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T, Q, I> ZSorted for  NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> NonOverlapping for NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> MOCProperties for NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> Iterator for NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{

  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    let new = if let Some(r) = self.it.next() {
      Some(Range {
        start: std::mem::replace(&mut self.start, r.end),
        end: r.start
      })
    } else if self.start == self.n_cells_max {
      None
    } else {
      Some(Range {
        start: std::mem::replace(&mut self.start, self.n_cells_max),
        end: self.n_cells_max,
      })
    };
    std::mem::replace(&mut self.curr, new)
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let cur = if self.curr.is_some() { 0 } else { 1 };
    let size_hint = self.it.size_hint();
    if let Some(n) = size_hint.1 {
      // +1 for cases in which the last range does not includes the upper bound
      (size_hint.0 + cur, Some(n  + cur + 1))
    } else {
      (size_hint.0 + cur, None)
    }
  }
}

impl<T, Q, I> RangeMOCIterator<T> for NotRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    // self.it.peek_last()
    // if not upper bound, return range.end..upper_bound
    // else return none
    None
  }

}
