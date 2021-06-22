
use std::io::Write;
use std::slice;
use std::ops::Range;
use std::vec::{IntoIter};
use std::marker::PhantomData;

use crate::mocrange::MocRange;
use crate::ranges::Idx;
use crate::mocqty::MocQty;
use crate::mocranges::MocRanges;
use crate::mocell::{Cell, CellOrCellRange, CellRange};
use crate::deser;
use crate::mocells::{MocCellOrCellRanges, MocCells, CellOrCellRanges};

/// Returns the maximum depth of an item the implementor contains.
pub trait HasMaxDepth {
  fn depth_max(&self) -> u8;
}

/// Marker trait telling that ranges lower bound or cells or indices in the implementor are sorted
/// according to the natural Z-order curve order (no hierarchical order => cells of various depth
/// are mixed).
pub trait ZSorted {}

/// Marker trait telling that something the implementor contains non-overlapping items.
pub trait NonOverlapping {}

pub trait MOCProperties: HasMaxDepth + ZSorted + NonOverlapping { }

////////////////////////////
// Encours d'elobaration
/*pub trait MOC: MOCProperties {
  type Idx: Idx;
  type Qty: MocQty<T>;
  type RangeIt: Iterator<Item=Range<Self::Idx>>;
  /// If we need to reuse the MOC, two solutions: clone or implement MOC on &Self.
  /// Then use Cow<MOC> in operations (if needed).
  fn into_iter(self) -> RangeIt;
}
pub trait MOCIterator<T: Idx>: Iterator<Item=Range<Idx>> {
  type Qty: MocQty<T>;
}

struct SweepLine {
  left_it
  right_it
}
Iterator<Item=Some(enum(Left(Ragne), Right(Range), Common(Range)))>.. ?*/
/////////////////////////////

// HierarchMOC
// RangesMOC
// CellsMOC



/*pub trait CellMOCIterator<T, Q>: Sized + MOCProperties + Iterator<Item=MocCell<T, Q>>
  where
    T: Idx,
    Q: MocQty<T>,
{

}
impl<T, Q, R> CellMOCIterator<T, Q> for R
  where
    T: Idx,
    Q: MocQty<T>,
    R: Sized + MOCProperties + Iterator<Item=MocCell<T, Q>> {}
*/

pub trait CellMOCIterator<T: Idx>: Sized + MOCProperties + Iterator<Item=Cell<T>> {
  type Qty: MocQty<T>;
  fn cellranges(self) -> CellOrCellRangeMOCIteratorFromCells<T, Self::Qty, Self> {
    CellOrCellRangeMOCIteratorFromCells::new(self)
  }
  fn ranges(self) -> RangeMOCIteratorFromCells<T, Self::Qty, Self> {
    RangeMOCIteratorFromCells::new(self)
  }
  fn to_json_aladin<W: Write>(self, fold: Option<usize>, writer: W) -> std::io::Result<()> {
    deser::json::to_json_aladin(self, &fold, "", writer)
  }
}
pub trait CellMOCIntoIterator<T: Idx>: Sized {
  type Qty: MocQty<T>;
  type IntoCellMOCIter: CellMOCIterator<T, Qty=Self::Qty>;

  fn into_cell_moc_iter(self) -> Self::IntoCellMOCIter;
}



pub trait CellOrCellRangeMOCIterator<T: Idx>: Sized + MOCProperties + Iterator<Item=CellOrCellRange<T>> {
  type Qty: MocQty<T>;

  /// # WARNING
  /// - `use_offset=true` is not compatible with the current IVOA standard!
  fn to_ascii_ivoa<W: Write>(self, fold: Option<usize>, use_offset: bool, writer: W) -> std::io::Result<()> {
    deser::ascii::to_ascii_ivoa(self, &fold, use_offset, writer)
  }
  /// # WARNING
  /// - this is not compatible with the current IVOA standard!
  fn to_ascii_stream<W: Write>(self, use_offset: bool, writer: W) -> std::io::Result<()> {
    deser::ascii::to_ascii_stream(self, use_offset, writer)
  }

  fn ranges(self) -> RangeMOCIteratorFromCellOrCellRanges<T, Self::Qty, Self> {
    RangeMOCIteratorFromCellOrCellRanges::new(self)
  }
}
pub trait CellOrCellRangeMOCIntoIterator<T: Idx>: Sized {
  type Qty: MocQty<T>;
  type IntoCellOrCellRangeMOCIter: CellOrCellRangeMOCIterator<T, Qty=Self::Qty>;

  fn into_cellcellrange_moc_iter(self) -> Self::IntoCellOrCellRangeMOCIter;
}



// Convert Cell --> Ranges
// TODO: merge (Cell --> Ranges) and (CellOrCellRange --> Ranges) in a single obj using Into<Range> ?

/// Transforms a `CellMOCIterator` into a `RangeMOCIterator`.
pub struct RangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  it: R,
  curr: Option<Range<T>>,
}
impl<T, Q, R> RangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  fn new(mut it: R) -> RangeMOCIteratorFromCells<T, Q, R> {
    let curr: Option<Range<T>> = it.next().map(|e| MocRange::<T, Q>::from(e).0);
    RangeMOCIteratorFromCells {
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
      Some(Range { start: lstart, end: lend }) => {
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
}



// Convert CellOrCellRanges --> Ranges

/// Transforms a `CellOrCellRangeMOCIterator` into a `RangeMOCIterator`.
pub struct RangeMOCIteratorFromCellOrCellRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellOrCellRangeMOCIterator<T, Qty=Q>
{
  it: R,
  curr: Option<Range<T>>,
}
impl<T, Q, R> RangeMOCIteratorFromCellOrCellRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellOrCellRangeMOCIterator<T, Qty=Q>
{
  fn new(mut it: R) -> RangeMOCIteratorFromCellOrCellRanges<T, Q, R> {
    let curr: Option<Range<T>> = it.next().map(|e| MocRange::<T, Q>::from(e).0);
    RangeMOCIteratorFromCellOrCellRanges {
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
      Some(Range { start: lstart, end: lend }) => {
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
}



// Convert Cells --> CellOrCellRanges

/// Transforms a `CellMOCIterator` into a `CellOrCellRangeMOCIterator`.
pub struct CellOrCellRangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  it: R,
  curr: Option<Cell<T>>,
}

impl<T, Q, R> CellOrCellRangeMOCIteratorFromCells<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: CellMOCIterator<T, Qty=Q>
{
  fn new(mut it: R) -> CellOrCellRangeMOCIteratorFromCells<T, Q, R> {
    let curr = it.next();
    CellOrCellRangeMOCIteratorFromCells {
      it,
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
}

// RangeMOCIterator<T, Q, B>: Sized + MOCProperties + Iterator<Item=B>
// with B: Borrow<MocRange<T, Q>>
// Don't do this, we will clone as iterating if necessary

/*
/// Here we define a MOC as a simple sequence of ranges having the following properties:
/// * the range bounds are expressed at Q::MAX_DEPTH
/// * the list contains mutually exclusive ranges
/// * the list is sorted from the smaller to the largest bounds
pub trait RangeMOCIterator<T, Q>: Sized + MOCProperties + Iterator<Item=MocRange<T, Q>>
  where
    T: Idx,
    Q: MocQty<T>,
{
  /// Transforms this iterator over ranges into a struct implementing the trait `MOCIterator`.
  fn to_cells_iter(self) -> CellMOCIteratorFromRanges<T, Q, Self> {
    CellMOCIteratorFromRanges::new(self)
  }
}
impl<T, Q, R> RangeMOCIterator<T, Q> for R
  where
    T: Idx,
    Q: MocQty<T>,
    R: Sized + MOCProperties + Iterator<Item=MocRange<T, Q>> {}
*/
pub trait RangeMOCIterator<T: Idx>: Sized + MOCProperties + Iterator<Item=Range<T>> {
  type Qty: MocQty<T>;

  fn cells(self) -> CellMOCIteratorFromRanges<T, Self::Qty, Self> {
    CellMOCIteratorFromRanges::new(self)
  }

}

pub trait RangeMOCIntoIterator<T: Idx>: Sized {
  type Qty: MocQty<T>;
  type IntoRangeMOCIter: RangeMOCIterator<T, Qty=Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter;
}

/// Transforms a `RangeMOCIterator` into a `CellMOCIterator`.
pub struct CellMOCIteratorFromRanges<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  it: R,
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

  fn new(mut it: R) -> CellMOCIteratorFromRanges<T, Q, R> {
    let curr = it.next().map(|range| range.into());
    let shift_dd = Q::shift_from_depth_max (it.depth_max()) as usize;
    let range_len_min = T::one() << shift_dd;
    let mut mask: T = From::from(Q::LEVEL_MASK);
    mask = mask.unsigned_shl(shift_dd as u32);
    CellMOCIteratorFromRanges {
      it,
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




// CELL/CELLRANGE MOC IMPLEMENTATION

/// A MOC made of a mix of (ordered and non-overlaping) cells and cells range
/// (a cell range is a range at the cell depth, while regular ranges are at the
/// largest possible depth and depends on type  `T`).
/// This is used for ASCII serialization/deserialization.
pub struct CellOrCellRangeMOC<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  ranges: MocCellOrCellRanges<T, Q>
}
impl<T: Idx, Q: MocQty<T>> CellOrCellRangeMOC<T, Q> {
  pub fn new(depth_max: u8, ranges: MocCellOrCellRanges<T, Q>) -> Self {
    Self { depth_max, ranges }
  }
  pub fn into_cellcellrange_moc_iter(self) -> CellOrCellRangeMocIter<T, Q> {
    CellOrCellRangeMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.into_iter(),
      _qty: PhantomData
    }
  }
  pub fn elems(self) -> CellOrCellRanges<T> {
    self.ranges.0
  }
  pub fn moc_elems(self) -> MocCellOrCellRanges<T, Q> {
    self.ranges
  }
}

impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellOrCellRangeMOC<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellOrCellRangeMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellOrCellRangeMOC<T, Q> { }

// - iterator
pub struct CellOrCellRangeMocIter<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: IntoIter<CellOrCellRange<T>>,
  _qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellOrCellRangeMocIter<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellOrCellRangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellOrCellRangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for CellOrCellRangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> Iterator for CellOrCellRangeMocIter<T, Q> {
  type Item = CellOrCellRange<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIterator<T> for CellOrCellRangeMocIter<T, Q> {
  type Qty = Q;
}
impl<T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIntoIterator<T> for CellOrCellRangeMOC<T, Q> {
  type Qty = Q;
  type IntoCellOrCellRangeMOCIter = CellOrCellRangeMocIter<T, Self::Qty>;

  fn into_cellcellrange_moc_iter(self) -> Self::IntoCellOrCellRangeMOCIter {
    CellOrCellRangeMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.into_iter(),
      _qty: PhantomData
    }
  }
}

// - ref iterator
pub struct CellOrCellRangeRefMocIter<'a, T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: slice::Iter<'a, CellOrCellRange<T>>,
  _qty: PhantomData<Q>,
}
impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for CellOrCellRangeRefMocIter<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<'a, T: Idx, Q: MocQty<T>> ZSorted for CellOrCellRangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for CellOrCellRangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for CellOrCellRangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> Iterator for CellOrCellRangeRefMocIter<'a, T, Q> {
  type Item = CellOrCellRange<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(|e| e.clone())
  }
  // Declaring size_hint, a 'collect' can directly allocate the right number of elements
  fn size_hint(&self) -> (usize, Option<usize>) {
    self.iter.size_hint()
  }
}
impl<'a, T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIterator<T> for CellOrCellRangeRefMocIter<'a, T, Q> {
  type Qty = Q;
}
impl<'a, T: Idx, Q: MocQty<T>> CellOrCellRangeMOCIntoIterator<T> for &'a CellOrCellRangeMOC<T, Q> {
  type Qty = Q;
  type IntoCellOrCellRangeMOCIter = CellOrCellRangeRefMocIter<'a, T, Self::Qty>;

  fn into_cellcellrange_moc_iter(self) -> Self::IntoCellOrCellRangeMOCIter {
    CellOrCellRangeRefMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.iter(),
      _qty: PhantomData
    }
  }
}

// CELL MOC IMPLEMENTATION
pub struct CellMOC<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  cells: MocCells<T, Q>
}
impl<T: Idx, Q: MocQty<T>> CellMOC<T, Q> {
  pub fn new(depth_max: u8, cells: MocCells<T, Q>) -> Self {
    Self { depth_max, cells }
  }
  pub fn len(&self) -> usize {
    self.cells.cells().cells().len()
  }
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellMOC<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for CellMOC<T, Q> { }


pub struct CellMocIter<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: IntoIter<Cell<T>>,
  _qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for CellMocIter<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for CellMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for CellMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for CellMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> Iterator for CellMocIter<T, Q> {
  type Item = Cell<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
}
impl<T: Idx, Q: MocQty<T>> CellMOCIterator<T> for CellMocIter<T, Q> {
  type Qty = Q;
}
impl<T: Idx, Q: MocQty<T>> CellMOCIntoIterator<T> for CellMOC<T, Q> {
  type Qty = Q;
  type IntoCellMOCIter = CellMocIter<T, Self::Qty>;

  fn into_cell_moc_iter(self) -> Self::IntoCellMOCIter {
    CellMocIter {
      depth_max: self.depth_max,
      iter: self.cells.0.0.into_iter(),
      _qty: PhantomData
    }
  }
}

// + Iterator<Item=Cell<T>>
// TODO: implementer l'iterateur, ...

pub struct CellRefMocIter<'a, T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: slice::Iter<'a, Cell<T>>,
  _qty: PhantomData<Q>,
}
impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for CellRefMocIter<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<'a, T: Idx, Q: MocQty<T>> ZSorted for CellRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for CellRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for CellRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> Iterator for CellRefMocIter<'a, T, Q> {
  type Item = Cell<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(|e| e.clone())
  }
}
impl<'a, T: Idx, Q: MocQty<T>> CellMOCIterator<T> for CellRefMocIter<'a, T, Q> {
  type Qty = Q;
}
impl<'a, T: Idx, Q: MocQty<T>> CellMOCIntoIterator<T> for &'a CellMOC<T, Q> {
  type Qty = Q;
  type IntoCellMOCIter = CellRefMocIter<'a, T, Self::Qty>;

  fn into_cell_moc_iter(self) -> Self::IntoCellMOCIter {
    CellRefMocIter {
      depth_max: self.depth_max,
      iter: self.cells.0.0.iter(),
      _qty: PhantomData
    }
  }
}



// RANGE MOC IMPLEMENTATION

/// A MOC made of (ordered and non-overlaping) ranges
#[derive(Debug, Clone)]
pub struct RangeMOC<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  ranges: MocRanges<T, Q>
}
impl<T: Idx, Q: MocQty<T>> RangeMOC<T, Q> {
  pub fn new(depth_max: u8, ranges: MocRanges<T, Q>) -> Self {
    Self {depth_max, ranges }
  }
  /// Returns the number of ranges the MOC contains
  pub fn len(&self) -> usize {
    self.ranges.0.0.len()
  }
  pub fn moc_ranges(&self) -> &MocRanges<T, Q> {
    &self.ranges
  }
  pub fn into_moc_ranges(self) -> MocRanges<T, Q> {
    self.ranges
  }
  /*pub fn into_range_moc_iter(self) -> LazyRangeMOCIter<T, Q, IntoIter<Range<T>>> {
    LazyRangeMOCIter::new(self.depth_max, self.ranges.0.0.into_iter())
  }*/

  /*pub fn range_moc_iter(&self) -> LazyRangeMOCVecIter<'_, H> {
    LazyRangeMOCVecIter::new(self.depth_max, self.ranges.iter())
  }*/
  /*pub fn into_cells_iter(self) -> CellMOCIteratorFromRanges<T, Q, Self> {
    CellMOCIteratorFromRanges::new(self)
  }*/
  /*pub fn to_cells_iter(&self) -> CellMOCIteratorFromRanges<T, Q, Self> {
    CellMOCIteratorFromRanges::new(self)
  }*/
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for RangeMOC<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for RangeMOC<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for RangeMOC<T, Q> { }


pub struct RangeMocIter<T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: IntoIter<Range<T>>,
  _qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>> HasMaxDepth for RangeMocIter<T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>> ZSorted for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> NonOverlapping for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> MOCProperties for RangeMocIter<T, Q> { }
impl<T: Idx, Q: MocQty<T>> Iterator for RangeMocIter<T, Q> {
  type Item = Range<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next()
  }
  // Declaring size_hint, a 'collect' can directly allocate the right number of elements
  fn size_hint(&self) -> (usize, Option<usize>) {
    self.iter.size_hint()
  }
}
impl<T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for RangeMocIter<T, Q> {
  type Qty = Q;
}
impl<T: Idx, Q: MocQty<T>> RangeMOCIntoIterator<T> for RangeMOC<T, Q> {
  type Qty = Q;
  type IntoRangeMOCIter = RangeMocIter<T, Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter {
    RangeMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.0.0.into_iter(),
      _qty: PhantomData
    }
  }
}

pub struct RangeRefMocIter<'a, T: Idx, Q: MocQty<T>> {
  depth_max: u8,
  iter: slice::Iter<'a, Range<T>>,
  _qty: PhantomData<Q>,
}
impl<'a, T: Idx, Q: MocQty<T>> HasMaxDepth for RangeRefMocIter<'a, T, Q> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<'a, T: Idx, Q: MocQty<T>> ZSorted for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> NonOverlapping for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> MOCProperties for RangeRefMocIter<'a, T, Q> { }
impl<'a, T: Idx, Q: MocQty<T>> Iterator for RangeRefMocIter<'a, T, Q> {
  type Item = Range<T>;
  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(|e| e.clone())
  }
  // Declaring size_hint, a 'collect' can directly allocate the right number of elements
  fn size_hint(&self) -> (usize, Option<usize>) {
    self.iter.size_hint()
  }
}
impl<'a, T: Idx, Q: MocQty<T>> RangeMOCIterator<T> for RangeRefMocIter<'a, T, Q> {
  type Qty = Q;
}
impl<'a, T: Idx, Q: MocQty<T>> RangeMOCIntoIterator<T> for &'a RangeMOC<T, Q> {
  type Qty = Q;
  type IntoRangeMOCIter = RangeRefMocIter<'a, T, Self::Qty>;

  fn into_range_moc_iter(self) -> Self::IntoRangeMOCIter {
    RangeRefMocIter {
      depth_max: self.depth_max,
      iter: self.ranges.iter(),
      _qty: PhantomData
    }
  }
}

// NUniq MOC
pub struct NUniqMOC<T: Idx> {
  depth_max: u8,
  zsorted_nuniq: Vec<T>
}
impl<T: Idx> NUniqMOC<T> {
  fn new(depth_max: u8, zsorted_nuniq: Vec<T>) -> Self {
    Self { depth_max, zsorted_nuniq}
  }
  /*fn from_unsorted(depth_max: u8, mut unsorted_nuniq: Vec<T>) -> Self {
    unsorted_nuniq.sort_by()
    Slef::new(depth_max, unsorted_nuniq)
  }*/
}





/*
/// A very basic, unchecked, Range MOC iterator.
/// We call it lazy because if the stored iterator is the result of a operation,
/// the operation will be performed only as we iterate over the LazyIterator: this
/// is an important point to be aware of when performing benches.
/// For an structure actually performing the operation, see `OwnedRangeMOC`.
pub struct LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  depth_max: u8,
  it: R,
  _t_t: PhantomData<T>,
  _t_q: PhantomData<Q>
}

impl<T, Q, R> LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  pub fn new(depth_max: u8, it: R) -> LazyRangeMOCIter<T, Q, R> {
    LazyRangeMOCIter { depth_max, it, _t_t: PhantomData, _t_q: PhantomData }
  }
}

impl<T, Q, R> HasMaxDepth for LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<T, Q, R> Iterator for LazyRangeMOCIter<T, Q, R>
  where
    T: Idx,
    Q: MocQty<T>,
    R: RangeMOCIterator<T, Qty=Q>
{
  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next()
  }
}
*/




/*
pub struct LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {
  depth_max: u8,
  iter: Iter<'a, HpxRange<H>>
}

impl<'a, H> LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {
  pub fn new(depth_max: u8, iter: Iter<'a, HpxRange<H>>) -> LazyRangeMOCVecIter<'a, H> {
    LazyRangeMOCVecIter{ depth_max, iter}
  }
}

impl<'a, H> HasMaxDepth for LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<'a, H> Iterator for LazyRangeMOCVecIter<'a, H>
  where H: HpxHash {

  type Item = HpxRange<H>;

  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(| HpxRange{ from, to} | HpxRange{ from: *from, to: *to})
  }
}
*/

#[cfg(test)]
mod tests {
  use crate::moc::{
    RangeMOC, RangeMOCIterator, RangeMOCIntoIterator,
    CellMOCIterator, CellOrCellRangeMOCIterator
  };
  use crate::mocranges::MocRanges;
  use crate::mocqty::Hpx;
  use crate::mocell::{Cell, CellRange, CellOrCellRange};
  use crate::hpxranges::HpxUniq2DepthIdxIter;
  use std::cmp::Ordering;

  #[test]
  fn test_range2cells_1() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![0..5])
    };
    let rit = (&rm).into_range_moc_iter();
    let v1: Vec<Cell<u64>> = rit.cells().collect();
    assert_eq!(v1, vec![
      Cell::new(28, 0),
      Cell::new(29, 4)
    ]);

    let v2: Vec<(i8, u64)> = HpxUniq2DepthIdxIter::new(rm.ranges).collect();
    assert_eq!(v2, vec![
      (28, 0),
      (29, 4)
    ]);
  }

  #[test]
  fn test_range2cells_2() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![2..8])
    };
    let rit = (&rm).into_range_moc_iter();
    let v1: Vec<Cell<u64>> = rit.cells().collect();
    assert_eq!(v1, vec![
      Cell::new(29, 2),
      Cell::new(29, 3),
      Cell::new(28, 1)
    ]);

    let v2: Vec<(i8, u64)> = HpxUniq2DepthIdxIter::new(rm.ranges).collect();
    assert_eq!(v2, vec![
      (28, 1),
      (29, 2),
      (29, 3)
    ]);
  }

  #[test]
  fn test_range2cells_3() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![
        0..5,
        6..59,
        78..6953,
        12458..55587,
        55787..65587
      ])
    };
    let rit = (&rm).into_range_moc_iter();
    let mut v1: Vec<Cell<u64>> = rit.cells().collect();
    // println!("{:?}", v1);
    v1.sort_by(
      |a, b| match a.depth.cmp(&b.depth) {
        Ordering::Less => Ordering::Less,
        Ordering::Greater => Ordering::Greater,
        Ordering::Equal => a.idx.cmp(&b.idx)
      });

    let v2: Vec<(i8, u64)> = HpxUniq2DepthIdxIter::new(rm.ranges).collect();
    // println!("{:?}", v2);
    assert_eq!(v1.len(), v2.len());
    for (Cell{depth, idx}, (depth2, idx2)) in v1.into_iter().zip(v2.into_iter())  {
      assert_eq!(depth, depth2 as u8);
      assert_eq!(idx, idx2);
    }
  }

  #[test]
  fn test_range2cellrange() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![2..8])
    };
    let res: Vec<CellOrCellRange<u64>> = (&rm).into_range_moc_iter().cells().cellranges().collect();
    assert_eq!(res, vec![
      CellOrCellRange::CellRange(CellRange::new(29, 2, 4)),
      CellOrCellRange::Cell(Cell::new(28, 1)),
    ]);

  }

  #[test]
  fn test_to_ascii() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![
        0..5,
        6..59,
        78..6953,
        12458..55587,
        55787..65587
      ])
    };
    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(Some(80), false, &mut sink)
      .unwrap();
    // println!("{}\n", &String::from_utf8_lossy(&sink));

    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(Some(80), true, &mut sink)
      .unwrap();
    //  println!("{}\n", &String::from_utf8_lossy(&sink));

    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(None, false, &mut sink)
      .unwrap();
    println!("{}\n", &String::from_utf8_lossy(&sink));

    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(None, true, &mut sink)
      .unwrap();
    //  println!("{}\n", &String::from_utf8_lossy(&sink));

    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_stream(false, &mut sink)
      .unwrap();
    //  println!("{}\n", &String::from_utf8_lossy(&sink));

    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_stream(true, &mut sink)
      .unwrap();
    //  println!("{}\n", &String::from_utf8_lossy(&sink));
  }

  #[test]
  fn test_to_json() {
    let rm = RangeMOC {
      depth_max: 29,
      ranges: MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![
        0..5,
        6..59,
        78..6953,
        12458..55587,
        55787..65587
      ])
    };
    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .to_json_aladin(Some(40), &mut sink)
      .unwrap();
    //  println!("{}\n", &String::from_utf8_lossy(&sink));

    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .to_json_aladin(None, &mut sink)
      .unwrap();
    //  println!("{}\n", &String::from_utf8_lossy(&sink));
  }

}