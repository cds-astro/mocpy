use crate::ranges::Idx;
use crate::mocqty::MocQty;
use crate::mocell::{Cell, CellOrCellRange};
use std::marker::PhantomData;

#[derive(Clone, Debug)]
pub struct Cells<T: Idx>(pub Vec<Cell<T>>);

#[derive(Debug)]
pub struct MocCells<T: Idx, Q: MocQty<T>>(pub Cells<T>, PhantomData<Q>);
impl<T: Idx, Q: MocQty<T>> MocCells<T, Q> {
  pub fn new(cells: Cells<T>) -> Self {
    Self(cells, PhantomData)
  }
}
// We could have chosen (Vec<MocCell<T, Q>)

#[derive(Debug)]
pub struct CellOrCellRanges<T: Idx>(pub Vec<CellOrCellRange<T>>);

#[derive(Debug)]
pub struct MocCellOrCellRanges<T: Idx, Q: MocQty<T>>(pub CellOrCellRanges<T>, PhantomData<Q>);
impl<T: Idx, Q: MocQty<T>> MocCellOrCellRanges<T, Q> {
  pub fn new(cells_or_cellranges: CellOrCellRanges<T>) -> Self {
    Self(cells_or_cellranges, PhantomData)
  }
}