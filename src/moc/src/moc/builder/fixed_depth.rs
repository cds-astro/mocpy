//! Builder in which we add cells at a given depth.

use std::slice;
use std::ops::Range;
use std::vec::IntoIter;
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::elemset::range::MocRanges;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator, range::RangeMOC, RangeMOCIntoIterator};
use crate::moc::range::op::or::or;

pub struct FixedDepthMocBuilder<T: Idx, Q: MocQty<T>> {
  depth: u8,
  buff: Vec<T>,
  sorted: bool,
  moc: Option<RangeMOC<T, Q>>,
}

impl<T: Idx, Q: MocQty<T>> FixedDepthMocBuilder<T, Q> {

  pub fn new(depth: u8, buf_capacity: Option<usize>) -> Self {
    FixedDepthMocBuilder {
      depth,
      buff: Vec::with_capacity(buf_capacity.unwrap_or(100_000)),
      sorted: true,
      moc: None
    }
  }

  pub fn from(buf_capacity: Option<usize>, moc: RangeMOC<T, Q>) -> Self {
    FixedDepthMocBuilder {
      depth: moc.depth_max(),
      buff: Vec::with_capacity(buf_capacity.unwrap_or(100_000)),
      sorted: true,
      moc: Some(moc)
    }
  }

  pub fn into_moc(mut self) -> RangeMOC<T, Q> {
    (&mut self).drain_buffer();
    let depth = self.depth;
    self.moc.unwrap_or_else(|| RangeMOC::new(depth, Default::default()))
  }

  pub fn into_moc_v2(mut self) -> RangeMOC<T, Q> {
    (&mut self).drain_buffer_v2();
    let depth = self.depth;
    self.moc.unwrap_or_else(|| RangeMOC::new(depth, Default::default()))
  }


  pub fn push(&mut self, idx: T) {
    if let Some(h) = self.buff.last() {
      if *h == idx {
        return;
      } else if self.sorted && *h > idx {
        self.sorted = false;
      }
    }
    self.buff.push(idx);
    if self.buff.len() == self.buff.capacity() {
      self.drain_buffer();
    }
  }

  // Bench against push()
  pub fn push_v2(&mut self, idx: T) {
    if let Some(h) = self.buff.last() {
      if *h == idx {
        return;
      } else if self.sorted && *h > idx {
        self.sorted = false;
      }
    }
    self.buff.push(idx);
    if self.buff.len() == self.buff.capacity() {
      self.drain_buffer_v2();
    }
  }


  /* do not use the internal buffer
  pub fn push_all(&mut self, idx: &[T]) {
    // sort
  }*/

  fn drain_buffer(&mut self) {
    if !self.sorted {
      // Sort without removing duplicates
      self.buff.sort_unstable();
    }
    let new_moc = self.buff_to_moc();
    self.clear_buff();
    let merged_moc = if let Some(prev_moc) = &self.moc {
      prev_moc.or(&new_moc)
    } else {
      new_moc
    };
    self.moc.replace(merged_moc);
  }

  fn drain_buffer_v2(&mut self) {
    if !self.sorted {
      // Sort without removing duplicates
      self.buff.sort_unstable();
    }
    let it = OrderedFixedDepthCellsToRanges::new(self.depth, self.buff.iter());
    let merged_moc = if let Some(prev_moc) = &self.moc {
      let or = or(prev_moc.into_range_moc_iter(),it);
      RangeMOC::new(self.depth, or.collect())
    } else {
      RangeMOC::new(self.depth, it.collect())
    };
    self.clear_buff();
    self.moc.replace(merged_moc);
  }

  // We make a copy of (a part of) the buffer
  fn buff_to_moc(&self) -> RangeMOC<T, Q> {
    let shift = Q::shift_from_depth_max(self.depth) as u32;
    let mut ranges = Vec::with_capacity(self.buff.len());
    // We assume here that the buffer is ordered, but may contains duplicates
    let mut it = self.buff.iter();
    if let Some(from) = it.next() {
      let mut from = *from;
      let mut to = from + T::one();
      for curr in it {
        if *curr == to {
          to += T::one();
        } else if *curr > to {
          ranges.push(from.unsigned_shl(shift)..to.unsigned_shl(shift));
          from = *curr;
          to = *curr + T::one();
        } else {
          debug_assert_eq!(*curr, to - T::one());
        }
      }
      ranges.push(from.unsigned_shl(shift)..to.unsigned_shl(shift));
    }
    RangeMOC::new(self.depth, MocRanges::new_unchecked(ranges))
  }

  fn clear_buff(&mut self) {
    self.sorted = true;
    self.buff.clear();
  }
}

/// The purpose of this struct is to avoid a copy, even if operations on
/// iterator is slower than operations on owned types.
pub struct OrderedFixedDepthCellsToRanges<'a, T: Idx, Q: MocQty<T>> {
  /// The max depth, needed to implement RangeMOCIterator
  depth: u8,
  /// The len of a range a cell represents
  shift: u32,
  /// Iterator over a sorted
  iter: slice::Iter<'a, T>,
  /// Current starting idx (inclusive)
  from: T,
  /// Current ending idx (exclusive)
  to: T,
  /// Tell that all elements have akready been returned
  depleted: bool,
  /// Q type
  _q_type: PhantomData<Q>
}

impl<'a, T: Idx, Q: MocQty<T>> OrderedFixedDepthCellsToRanges<'a, T, Q> {
  pub fn new(depth: u8, mut iter: slice::Iter<'a, T>) -> Self {
    let shift = Q::shift_from_depth_max(depth) as u32;
    let (from, to, depleted) = if let Some(v) = iter.next() {
      (v.unsigned_shl(shift), (*v + T::one()).unsigned_shl(shift), false)
    } else {
      (T::zero(), T::one(), true)
    };
    OrderedFixedDepthCellsToRanges {
      depth,
      shift,
      iter,
      from,
      to,
      depleted,
      _q_type: PhantomData
    }
  }
}

impl<'a, T: Idx, Q: MocQty<T>> Iterator for OrderedFixedDepthCellsToRanges<'a, T, Q> {
  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    while let Some(curr) = self.iter.next() {
      if *curr == self.to {
        self.to += T::one();
      } else if *curr > self.to {
        let range = self.from.unsigned_shl(self.shift)..self.to.unsigned_shl(self.shift);
        self.from = *curr;
        self.to = *curr + T::one();
        return Some(range)
      } else {
        debug_assert_eq!(*curr, self.to - T::one());
      }
    }
    if !self.depleted {
      self.depleted = true;
      let range = self.from.unsigned_shl(self.shift)..self.to.unsigned_shl(self.shift);
      Some(range)
    } else {
      None
    }
  }
}

impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for OrderedFixedDepthCellsToRanges<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth
  }
}
impl<'a, T: Idx, Q: MocQty<T>> ZSorted for OrderedFixedDepthCellsToRanges<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for OrderedFixedDepthCellsToRanges<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for OrderedFixedDepthCellsToRanges<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for OrderedFixedDepthCellsToRanges<'a, T, Q> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    None
  }
}




/// The purpose of this struct is to avoid a copy, even if operations on
/// iterator is slower than operations on owned types.
pub struct OwnedOrderedFixedDepthCellsToRanges<T: Idx, Q: MocQty<T>> {
  /// The max depth, needed to implement RangeMOCIterator
  depth: u8,
  /// The len of a range a cell represents
  shift: u32,
  /// Iterator over a sorted
  iter: IntoIter<T>,
  /// Current starting idx (inclusive)
  from: T,
  /// Current ending idx (exclusive)
  to: T,
  /// Tell that all elements have akready been returned
  depleted: bool,
  /// Q type
  _q_type: PhantomData<Q>
}

impl<T: Idx, Q: MocQty<T>> OwnedOrderedFixedDepthCellsToRanges<T, Q> {
  pub fn new(depth: u8, mut iter: IntoIter<T>) -> Self {
    let shift = Q::shift_from_depth_max(depth) as u32;
    let (from, to, depleted) = if let Some(v) = iter.next() {
      (v.unsigned_shl(shift), (v + T::one()).unsigned_shl(shift), false)
    } else {
      (T::zero(), T::one(), true)
    };
    OwnedOrderedFixedDepthCellsToRanges {
      depth,
      shift,
      iter,
      from,
      to,
      depleted,
      _q_type: PhantomData
    }
  }
}

impl<T: Idx, Q: MocQty<T>> Iterator for OwnedOrderedFixedDepthCellsToRanges<T, Q> {
  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    while let Some(curr) = self.iter.next() {
      if curr == self.to {
        self.to += T::one();
      } else if curr > self.to {
        let range = self.from.unsigned_shl(self.shift)..self.to.unsigned_shl(self.shift);
        self.from = curr;
        self.to = curr + T::one();
        return Some(range)
      } else {
        debug_assert_eq!(curr, self.to - T::one());
      }
    }
    if !self.depleted {
      self.depleted = true;
      let range = self.from.unsigned_shl(self.shift)..self.to.unsigned_shl(self.shift);
      Some(range)
    } else {
      None
    }
  }
}

impl<T: Idx, Q: MocQty<T>> HasMaxDepth for OwnedOrderedFixedDepthCellsToRanges<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for OwnedOrderedFixedDepthCellsToRanges<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for OwnedOrderedFixedDepthCellsToRanges<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for OwnedOrderedFixedDepthCellsToRanges<T, Q> { }
impl<T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for OwnedOrderedFixedDepthCellsToRanges<T, Q> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    None
  }
}