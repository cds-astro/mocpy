
use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};

/// Performs a `degradation` on the input iterator of ranges.
pub fn degrade<T, Q, I>(it: I, new_depth: u8) -> DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>
{
  DegradeRangeIter::new(it, new_depth)
}

/// Performs a `degradation` on the input iterator of ranges.
pub struct DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  new_depth: u8,
  one_at_new_depth: T,
  rm_bits_mask: T,
  bits_to_be_rm_mask: T,
  it: I,
  curr: Option<Range<T>>
}

impl<T, Q, I> DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{

  pub fn new(mut it: I, new_depth: u8) -> DegradeRangeIter<T, Q, I>  {
    if new_depth < it.depth_max() {
      let shift = Q::shift_from_depth_max(new_depth) as u32;
      let one_at_new_depth = T::one().unsigned_shl(shift);
      let rm_bits_mask = (!T::zero()).unsigned_shl(shift);
      let bits_to_be_rm_mask = !rm_bits_mask; // If (end & test_end_mask) == 0, do nothing, else + one
      let mut curr = it.next();
      if let Some(cr) = &mut curr {
        degrade_range(cr, one_at_new_depth, rm_bits_mask, bits_to_be_rm_mask)
      }
      DegradeRangeIter {
        new_depth,
        one_at_new_depth,
        rm_bits_mask,
        bits_to_be_rm_mask,
        it,
        curr,
      }
    } else {
      // Do nothing, just change the depth.
      let new_depth = it.depth_max();
      let curr = it.next();
      DegradeRangeIter {
        new_depth,
        one_at_new_depth: T::one(), 
        rm_bits_mask: !T::zero(),
        bits_to_be_rm_mask: T::zero(),
        it,
        curr
      }
    }
  }
}

#[inline(always)]
fn degrade_range<T: Idx>(range: &mut Range<T>, one_at_new_depth: T, rm_bits_mask: T, bits_to_be_rm_mask: T) {
  range.start &= rm_bits_mask;
  if range.end & bits_to_be_rm_mask != T::zero() {
    range.end = (range.end & rm_bits_mask) + one_at_new_depth;
  }
}

impl<T, Q, I> HasMaxDepth for  DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  fn depth_max(&self) -> u8 {
    self.new_depth
  }
}
impl<T, Q, I> ZSorted for  DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> NonOverlapping for DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> MOCProperties for DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{ }
impl<T, Q, I> Iterator for DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{

  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    if let Some(cr) = &mut self.curr {
      let mut next = self.it.next();
      while let Some(nr) = &mut next {
        degrade_range(nr, self.one_at_new_depth, self.rm_bits_mask, self.bits_to_be_rm_mask);
        if nr.start > cr.end {
          break;
        } else {
          cr.end = nr.end;
          next = self.it.next();
        }
      }
      std::mem::replace(&mut self.curr, next)
    } else {
      None
    }
  }

  /*fn size_hint(&self) -> (usize, Option<usize>) {
  }*/
}

impl<T, Q, I> RangeMOCIterator<T> for DegradeRangeIter<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    // self.it.peek_last()
    // TODO: peek last, degrade and store the result to provide it here
    None
  }

}
