
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::elem::cell::Cell;

#[derive(Clone, Debug)]
pub struct Cells<T: Idx>(pub Box<[Cell<T>]>);
impl<T: Idx> Cells<T> {
  pub fn new(cells: Vec<Cell<T>>) -> Self {
    Self(cells.into_boxed_slice())
  }
  pub fn cells(&self) -> &[Cell<T>] { &self.0 }
}

#[derive(Debug)]
pub struct MocCells<T: Idx, Q: MocQty<T>>(pub Cells<T>, PhantomData<Q>);
impl<T: Idx, Q: MocQty<T>> MocCells<T, Q> {
  pub fn new(cells: Cells<T>) -> Self {
    Self(cells, PhantomData)
  }
  pub fn cells(&self) -> &[Cell<T>] { self.0.cells() }
}
// We could have chosen (Vec<MocCell<T, Q>)

