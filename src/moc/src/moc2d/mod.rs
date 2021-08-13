
use std::io::Write;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{
  ZSorted, NonOverlapping, 
  CellOrCellRangeMOCIterator, CellMOCIterator, RangeMOCIterator
};

use crate::deser::ascii::{AsciiError, moc2d_to_ascii_ivoa};
use crate::deser::json::cellmoc2d_to_json_aladin;
use crate::moc2d::range::{RangeMOC2, RangeMOC2Elem};

pub mod cell;
pub mod range;
pub mod cellcellrange;
pub mod adapters;

/// Returns the maximum depth of an item the implementor contains.
pub trait HasTwoMaxDepth {
  fn depth_max_1(&self) -> u8;
  fn depth_max_2(&self) -> u8;
}
/// Must have all good properties
pub trait MOC2Properties: HasTwoMaxDepth + ZSorted + NonOverlapping { }


// Traits related to CellOrCellRange (ASCII format)

pub trait CellOrCellRangeMOC2ElemIt<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  type It1: CellOrCellRangeMOCIterator<T, Qty=Q>;
  type It2: CellOrCellRangeMOCIterator<U, Qty=R>;
  fn cellcellrange_mocs_it(self) -> (Self::It1, Self::It2);
}

pub trait CellOrCellRangeMOC2Iterator<
  T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
  K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>
>: Sized
  + MOC2Properties
  + Iterator<Item=K>
{
  /// # WARNING
  /// - `use_offset=true` is not compatible with the current IVOA standard!
  fn to_ascii_ivoa<W: Write>(self, fold: Option<usize>, use_offset: bool, writer: W) -> Result<(), AsciiError> {
    moc2d_to_ascii_ivoa(self, &fold, use_offset, writer)
  }
}

pub trait CellOrCellRangeMOC2IntoIterator<
  T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
  K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>
>: Sized {
  type IntoCellOrCellRangeMOC2Iter: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K>;

  fn into_cellcellrange_moc2_iter(self) -> Self::IntoCellOrCellRangeMOC2Iter;
}

// Traits related to Cell (JSON format)

pub trait CellMOC2ElemIt<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  type It1: CellMOCIterator<T, Qty=Q>;
  type It2: CellMOCIterator<U, Qty=R>;
  fn cell_mocs_it(self) -> (Self::It1, Self::It2);
}

pub trait CellMOC2Iterator<
  T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
  K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>
>: Sized
  + MOC2Properties
  + Iterator<Item=K>
{
  /*/// # WARNING
  /// - `use_offset=true` is not compatible with the current IVOA standard!
  fn to_ascii_ivoa<W: Write>(self, fold: Option<usize>, use_offset: bool, writer: W) -> Result<(), AsciiError> {
    deser::ascii::moc2d_to_ascii_ivoa(self, &fold, use_offset, writer)
  }*/

  fn to_json_aladin<W: Write>(self, fold: &Option<usize>, writer: W) -> std::io::Result<()> {
    cellmoc2d_to_json_aladin(self, fold, writer)
  }
}

pub trait CellMOC2IntoIterator<
  T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
  K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>
>: Sized {
  type IntoCellMOC2Iter: CellMOC2Iterator<T, Q, I, U, R, J, K>;

  fn into_cell_moc2_iter(self) -> Self::IntoCellMOC2Iter;
}

// Traits related to Range (FITS format)

pub trait RangeMOC2ElemIt<T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>> {
  type It1: RangeMOCIterator<T, Qty=Q>;
  type It2: RangeMOCIterator<U, Qty=R>;

  fn range_mocs_it(self) -> (Self::It1, Self::It2);
}

pub trait RangeMOC2Iterator<
  T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
  K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>
>: Sized
  + MOC2Properties
  + Iterator<Item=K>
{
  fn into_range_moc2(self) -> RangeMOC2<T, Q, U, R> {
    let depth_max_l = self.depth_max_1();
    let depth_max_r = self.depth_max_2();
    let elems: Vec<RangeMOC2Elem<T, Q, U, R>> = self.map(|e| {
        let (it_l, it_r) = e.range_mocs_it();
        RangeMOC2Elem::new(it_l.into_range_moc(), it_r.into_range_moc())
      })
      .collect();
    RangeMOC2::new(depth_max_l, depth_max_r, elems)
  }
  /*
  fn to_fits_ivoa<W: Write>(
    self,
    moc_id: Option<String>,
    moc_type: Option<MocType>,
    mut writer: W
  ) -> Result<(), FitsError> {
    THE FOLLOWING METHOD IS NOT GENERICAL ENOUGH...
    ranges2d_to_fits_ivoa(self, moc_id, moc_type, writer)
  }*/
}

pub trait RangeMOC2IntoIterator<
  T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
  K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>
>: Sized {
  type IntoRangeMOC2Iter: RangeMOC2Iterator<T, Q, I, U, R, J, K>;

  fn into_range_moc2_iter(self) -> Self::IntoRangeMOC2Iter;
}
