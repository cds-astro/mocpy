
use std::ops::Range;
use std::vec::IntoIter;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};
use std::marker::PhantomData;

/// Build a range MOC from unordered, possibly overlapping ranges.
/// WARNING: no degrade operation performed, we trust the provided `depth_max`
/// without checking it!
pub fn merge_random<T, Q>(depth_max: u8, mut ranges: Vec<Range<T>>)
  -> MergeIterator<T, Q, IntoIter<Range<T>>>
  where
    T: Idx,
    Q: MocQty<T>
{
  ranges.sort_unstable_by(|a, b| a.start.cmp(&b.start));
  merge_sorted(depth_max, ranges.into_iter())
}

/// Decorates the given iterator with an iterator that panics (while iterating) if the input
/// iterator is not made of sorted, non overlaping ranges.
/// WARNING: no degrade operation performed, we trust the provided `depth_max`
/// without checking it!
pub fn merge_sorted<T, Q, I>(depth_max: u8, it: I) -> MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{
  // debug_assert is_sorted (so far nightly: https://doc.rust-lang.org/std/vec/struct.Vec.html#method.is_sorted_by)
  MergeIterator::new(depth_max, it)
}

/// Iterator decorator made to ensure that the decorated iterator returns sorted non-overlapping
/// ranges.
/// If it is not the case, a call to the `next` method `panics`!
pub struct MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{
  depth: u8,
  it: I,
  curr: Option<Range<T>>,
  _q_type: PhantomData<Q>,
}

impl<T, Q, I> MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{

  pub fn new(depth_max: u8, mut it: I) -> MergeIterator<T, Q, I>  {
    let curr = it.next();
    MergeIterator {
      depth: depth_max,
      it,
      curr,
      _q_type: PhantomData
    }
  }
}

impl<T, Q, I> HasMaxDepth for  MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{
  fn depth_max(&self) -> u8 {
    self.depth
  }
}
impl<T, Q, I> ZSorted for  MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{ }
impl<T, Q, I> NonOverlapping for MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{ }
impl<T, Q, I> MOCProperties for MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{ }
impl<T, Q, I> Iterator for MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{

  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    if let Some(prev) =  &mut self.curr {
      let mut curr_opt = self.it.next();
      while let Some(curr) = &curr_opt {
        debug_assert!(curr.start >= prev.start);
        if curr.start <= prev.end {
          prev.end = prev.end.max(curr.end);
          curr_opt = self.it.next();
        } else {
          break;
        }
      }
      std::mem::replace(&mut self.curr, curr_opt)
    } else {
      None
    }
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let mut sh = self.it.size_hint();
    sh.0 = if sh.0 > 0 { 1 } else { 0 };
    sh.1 = sh.1.map(|v| v + 1);
    sh
  }
}

impl<T, Q, I> RangeMOCIterator<T> for MergeIterator<T, Q, I>
  where
    T: Idx,
    Q: MocQty<T>,
    I: Iterator<Item=Range<T>>,
{
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    None
  }

}


#[cfg(test)]
mod tests {

  use core::ops::Range;

  use crate::qty::Hpx;
  use crate::moc::range::op::merge::merge_sorted;

  #[test]
  fn test_merge() {
    let input: Vec<Range<u32>> = vec![3..4, 6..10, 6..12, 8..14, 16..17];
    let res: Vec<Range<u32>> = merge_sorted::<u32, Hpx::<u32>, _>(0, input.into_iter()).collect();
    assert_eq!(
      res,
      vec![3..4, 6..14, 16..17]
    )
  }
}
