use crate::ranges::Idx;
use crate::mocqty::MocQty;
use std::marker::PhantomData;
use std::ops::Range;
use std::cmp::Ordering;

/// The motivation for this enum is the ASCII serialization which looks like:
/// > 3/3 10 4/16-18 22 5/19-20 17/222 28/123456789 29/
/// Mixing single cells and cells range.
/// This is usefull for Qty having a DIM > 1, because at DIM = 1 a cell is only divided in 2
/// (so we get a super cell instead of a range).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CellOrCellRange<T: Idx> {
  Cell(Cell<T>),
  CellRange(CellRange<T>),
}
impl<T: Idx> CellOrCellRange<T> {
  fn get_depth_idx_low(&self) -> (u8, &T) {
    match &self {
      CellOrCellRange::Cell(Cell{ depth, idx}) => (*depth, idx),
      CellOrCellRange::CellRange(CellRange{ depth, range}) => (*depth, &range.start),
    }
  }
  pub fn cmp<Q: MocQty<T>>(&self, other: &Self) -> Ordering {
    let (d1, i_l) = self.get_depth_idx_low();
    let (d2, i_r) = other.get_depth_idx_low();
    match d1.cmp(&d2) {
      Ordering::Equal => i_l.cmp(i_r),
      Ordering::Less => i_l.unsigned_shr(Q::shift(d2 - d1) as u32).cmp(i_r),
      Ordering::Greater => i_l.cmp(&i_r.unsigned_shr(Q::shift(d1 - d2) as u32)),
    }
  }
}


#[derive(Debug, Clone, PartialEq)]
pub enum MocCellOrCellRange<T: Idx, Q: MocQty<T>> {
  MocCell(MocCell<T, Q>),
  MocCellRange(MocCellRange<T, Q>),
}


/// A MOC cell, i.e. an index at a given depth.
/// Without attached quantities, we do not know the shift from one depth to another, and
/// so we cannot define a absolute order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cell<T: Idx> {
  pub depth: u8,
  pub idx: T,
}
impl<T: Idx> Cell<T> {
  pub fn new(depth: u8, idx: T) -> Self {
    Cell {depth, idx}
  }
  pub fn cmp<Q: MocQty<T>>(&self, other: &Self) -> Ordering {
    match self.depth.cmp(&other.depth) {
      Ordering::Equal => self.idx.cmp(&other.idx),
      Ordering::Less => self.idx.unsigned_shr(Q::shift(other.depth - self.depth) as u32).cmp(&other.idx),
      Ordering::Greater => self.idx.cmp(&other.idx.unsigned_shr(Q::shift(self.depth - other.depth) as u32)),
    }
  }
}

impl<T: Idx, Q: MocQty<T>> From<MocCell<T, Q>> for Cell<T> {
  fn from(cell: MocCell<T, Q>) -> Self {
    cell.0
  }
}

/// The order we define corresponds to the order of the lower bound of the cell at the deepest depth.
#[repr(transparent)] // To be able to transmute Vec<Cell<T>> into Vec<MocCell<T, _>>
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MocCell<T: Idx, Q: MocQty<T>>(pub Cell<T>, PhantomData<Q>);

impl<T: Idx, Q: MocQty<T>> Ord for MocCell<T, Q> {
  fn cmp(&self, other: &Self) -> Ordering {
    match self.0.depth.cmp(&other.0.depth) {
      Ordering::Equal => self.0.idx.cmp(&other.0.idx),
      Ordering::Less => self.0.idx.unsigned_shr(Q::shift(other.0.depth - self.0.depth) as u32).cmp(&other.0.idx),
      Ordering::Greater => self.0.idx.cmp(&other.0.idx.unsigned_shr(Q::shift(self.0.depth - other.0.depth) as u32)),
    }
  }
}
impl<T: Idx, Q: MocQty<T>> PartialOrd for MocCell<T, Q> {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

impl<T: Idx, Q: MocQty<T>> From<Cell<T>> for MocCell<T, Q> {
  fn from(cell: Cell<T>) -> Self {
    Self(cell, PhantomData)
  }
}


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
  pub fn cmp<Q: MocQty<T>>(&self, other: &Self) -> Ordering {
    match self.depth.cmp(&other.depth) {
      Ordering::Equal => self.range.start.cmp(&other.range.start),
      Ordering::Less => self.range.start.unsigned_shr(Q::shift(other.depth - self.depth) as u32).cmp(&other.range.start),
      Ordering::Greater => self.range.start.cmp(&other.range.start.unsigned_shr(Q::shift(self.depth - other.depth) as u32)),
    }
  }
}

/// The order we define corresponds to the order of the lower bound of the range at the deepest depth.
#[repr(transparent)] // To be able to transmute Vec<CellRange<T>> into Vec<MocCellRange<T, _>>
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MocCellRange<T: Idx, Q: MocQty<T>>(pub CellRange<T>, PhantomData<Q>);

impl<T: Idx, Q: MocQty<T>> Ord for MocCellRange<T, Q> {
  fn cmp(&self, other: &Self) -> Ordering {
    match self.0.depth.cmp(&other.0.depth) {
      Ordering::Equal => self.0.range.start.cmp(&other.0.range.start),
      Ordering::Less => self.0.range.start.unsigned_shr(Q::shift(other.0.depth - self.0.depth) as u32).cmp(&other.0.range.start),
      Ordering::Greater => self.0.range.start.cmp(&other.0.range.start.unsigned_shr(Q::shift(self.0.depth - other.0.depth) as u32)),
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