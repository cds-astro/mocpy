//! Builder in which we add ranges at the maximum depth.

use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, range::RangeMOC, RangeMOCIntoIterator};
use crate::moc::range::op::{
  or::or,
  merge::merge_sorted
};

pub struct RangeMocBuilder<T: Idx, Q: MocQty<T>> {
  depth: u8,
  buff: Vec<Range<T>>,
  sorted: bool,
  moc: Option<RangeMOC<T, Q>>,
}

impl<T: Idx, Q: MocQty<T>> RangeMocBuilder<T, Q> {

  pub fn new(depth: u8, buf_capacity: Option<usize>) -> Self {
    RangeMocBuilder {
      depth,
      buff: Vec::with_capacity(buf_capacity.unwrap_or(100_000)),
      sorted: true,
      moc: None
    }
  }

  pub fn from(buf_capacity: Option<usize>, moc: RangeMOC<T, Q>) -> Self {
    RangeMocBuilder {
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

  pub fn push(&mut self, range: Range<T>) {
    if let Some(Range { start, end }) = self.buff.last_mut() {
      if range.end < *start || *end < range.start {
        // both ranges do not overlap
        self.sorted &= range.end < *start;
        self.buff.push(range);
      } else {
        // merge overlaping ranges
        if range.start < *start {
          self.sorted = false; // we could try to look a previous ranges to merge them...
          *start = range.start;
        }
        *end = range.end.min(*end);
      }
    } else {
      self.buff.push(range);
    }
    if self.buff.len() == self.buff.capacity() {
      self.drain_buffer();
    }
  }

  fn drain_buffer(&mut self) {
    if !self.sorted {
      // Sort without removing duplicates
      self.buff.sort_unstable_by(|a, b| a.start.cmp(&b.start));
    }
    let new_moc_it = merge_sorted(self.depth, self.buff.drain(..));
    self.sorted = true;
    let merged_moc = if let Some(prev_moc) = &self.moc {
      let or = or(prev_moc.into_range_moc_iter(), new_moc_it);
      RangeMOC::new(self.depth, or.collect())
    } else {
      RangeMOC::new(self.depth, new_moc_it.collect())
    };
    self.moc.replace(merged_moc);
  }

  fn buff_to_moc(&mut self) -> RangeMOC<T, Q> {
    self.drain_buffer();
    self.moc.take().unwrap_or_else(|| RangeMOC::new(self.depth, Default::default()))
  }

}

