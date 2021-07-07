
use std::ops::Range;
use std::cmp::Ordering;
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::qty::MocQty;

/// A range of MOC cells, the range being expressed at the CellRange depth (and not at
/// the deepest possible depth like in a MocRange).
/// This is usefull for Qty having a DIM > 1, because at DIM = 1 a cell is only divided in 2
/// (so we get a super cell instead of a range).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CellRange<T: Idx> {
    pub depth: u8,
    pub range: Range<T>,
}
impl<T: Idx> CellRange<T> {
    pub fn new(depth: u8, start: T, end: T) -> Self {
        CellRange { depth, range: Range{ start, end}}
    }
    /// Comparison independent from the hierarchy, i.e. like a deepest level comparison.
    pub fn flat_cmp<Q: MocQty<T>>(&self, other: &Self) -> Ordering {
        match self.depth.cmp(&other.depth) {
            Ordering::Equal => self.range.start.cmp(&other.range.start),
            Ordering::Less => self.range.start.unsigned_shl(Q::shift(other.depth - self.depth) as u32).cmp(&other.range.start),
            Ordering::Greater => self.range.start.cmp(&other.range.start.unsigned_shl(Q::shift(self.depth - other.depth) as u32)),
        }
    }
}

/// The order we define corresponds to the order of the lower bound of the range at the deepest depth.
#[repr(transparent)] // To be able to transmute Vec<CellRange<T>> into Vec<MocCellRange<T, _>>
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MocCellRange<T: Idx, Q: MocQty<T>>(pub CellRange<T>, PhantomData<Q>);

impl<T: Idx, Q: MocQty<T>> Ord for MocCellRange<T, Q> {
    /// Comparison independent from the hierarchy, i.e. like a deepest level comparison.
    fn cmp(&self, other: &Self) -> Ordering {
        match self.0.depth.cmp(&other.0.depth) {
            Ordering::Equal => self.0.range.start.cmp(&other.0.range.start),
            Ordering::Less => self.0.range.start.unsigned_shl(Q::shift(other.0.depth - self.0.depth) as u32).cmp(&other.0.range.start),
            Ordering::Greater => self.0.range.start.cmp(&other.0.range.start.unsigned_shl(Q::shift(self.0.depth - other.0.depth) as u32)),
        }
    }
}
impl<T: Idx, Q: MocQty<T>> PartialOrd for MocCellRange<T, Q> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: Idx, Q: MocQty<T>> From<CellRange<T>> for MocCellRange<T, Q> {
    fn from(cell_range: CellRange<T>) -> Self {
        Self(cell_range, PhantomData)
    }
}