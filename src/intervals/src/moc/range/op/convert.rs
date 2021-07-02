
use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};
use std::marker::PhantomData;

/// Iterator decorator made to ensure that the decorated iterator returns sorted non-overlapping
/// ranges.
/// If it is not the case, a call to the `next` method `panics`!
pub struct ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{
  last: Option<Range<U>>,
  it: I,
  _t_type: PhantomData<T>,
  _q_type: PhantomData<Q>,
  _r_type: PhantomData<R>,
}

impl<T, Q, I, U, R> ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{

  pub fn new(it: I) -> ConvertIterator<T, Q, I, U, R>  {
    let last = it.peek_last().map(|r| r.start.convert()..r.end.convert());
    ConvertIterator {
      last,
      it,
      _t_type: PhantomData,
      _q_type: PhantomData,
      _r_type: PhantomData
    }
  }
}

impl<T, Q, I, U, R> HasMaxDepth for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T, Q, I, U, R> ZSorted for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{ }
impl<T, Q, I, U, R> NonOverlapping for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{ }
impl<T, Q, I, U, R> MOCProperties for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{ }
impl<T, Q, I, U, R> Iterator for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{

  type Item = Range<U>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next().map(|r| r.start.convert()..r.end.convert())
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let mut sh = self.it.size_hint();
    sh.0 += 1;
    sh.1 = sh.1.map(|v| v + 1);
    sh
  }
}

impl<T, Q, I, U, R> RangeMOCIterator<U> for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{
  type Qty = R;

  fn peek_last(&self) -> Option<&Range<U>> {
    self.last.as_ref()
  }
}