
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::elem::cellcellrange::CellOrCellRange;

#[derive(Debug)]
pub struct CellOrCellRanges<T: Idx>(pub Box<[CellOrCellRange<T>]>);
impl<T: Idx> CellOrCellRanges<T> {
    pub fn new(elems: Vec<CellOrCellRange<T>>) -> Self {
        Self(elems.into_boxed_slice())
    }
    pub fn elems(&self) -> &[CellOrCellRange<T>] { &self.0 }
}

#[derive(Debug)]
pub struct MocCellOrCellRanges<T: Idx, Q: MocQty<T>>(pub CellOrCellRanges<T>, PhantomData<Q>);
impl<T: Idx, Q: MocQty<T>> MocCellOrCellRanges<T, Q> {
    pub fn new(cells_or_cellranges: CellOrCellRanges<T>) -> Self {
        Self(cells_or_cellranges, PhantomData)
    }
    pub fn elems(&self) -> &[CellOrCellRange<T>] { self.0.elems() }
}