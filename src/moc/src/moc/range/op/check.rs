
use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};

/// Decorates the given iterator with an iterator that panics (while iterating) if the input
/// iterator is not made of sorted, non overlaping ranges.
pub fn check<T, Q, I>(it: I) -> CheckedIterator<T, Q, I>
  where
      T: Idx,
      Q: MocQty<T>,
      I: RangeMOCIterator<T, Qty=Q>
{
  CheckedIterator::new(it)
}

/// Iterator decorator made to ensure that the decorated iterator returns sorted non-overlapping
/// ranges.
/// If it is not the case, a call to the `next` method `panics`!
pub struct CheckedIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>, 
{
  it: I,
  curr: Option<Range<T>>,
}

impl<T, Q, I> CheckedIterator<T, Q, I> 
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>, 
{

  pub fn new(mut it: I) -> CheckedIterator<T, Q, I>  {
    let curr = it.next();
    CheckedIterator {
      it,
      curr,
    }
  }
}

impl<T, Q, I> HasMaxDepth for  CheckedIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T, Q, I> ZSorted for  CheckedIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> NonOverlapping for CheckedIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> MOCProperties for CheckedIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> Iterator for CheckedIterator<T, Q, I> 
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>, 
{

  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    let prev = std::mem::replace(&mut self.curr, self.it.next());
    if let (Some(l), Some(r)) = (&prev, &self.curr) {
      assert!(l.end < r.start) // not equals, else ranges should have been merged!
    }
    prev
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let mut sh = self.it.size_hint();
    sh.0 += 1;
    sh.1 = sh.1.map(|v| v + 1);
    sh
  }
}

impl<T, Q, I> RangeMOCIterator<T> for CheckedIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    self.it.peek_last()
  }

}
