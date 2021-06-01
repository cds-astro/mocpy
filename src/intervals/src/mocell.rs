use crate::ranges::Idx;
use crate::mocqty::MocQty;
use std::marker::PhantomData;
use std::ops::Range;

#[derive(Debug, Clone, PartialEq)]
pub struct Cell<T: Idx> {
  pub depth: u8,
  pub idx: T,
}

impl<T: Idx> Cell<T> {
  pub fn new(depth: u8, idx: T) -> Self {
    Cell {depth, idx}
  }
}

impl<T: Idx, Q: MocQty<T>> From<MocCell<T, Q>> for Cell<T> {
  fn from(cell: MocCell<T, Q>) -> Self {
    cell.0
  }
}



#[derive(Debug, Clone)]
pub struct MocCell<T: Idx, Q: MocQty<T>>(pub Cell<T>, PhantomData<Q>);

impl<T: Idx, Q: MocQty<T>> From<Cell<T>> for MocCell<T, Q> {
  fn from(cell: Cell<T>) -> Self {
    Self(cell, PhantomData)
  }
}

/*impl<T: Idx, Q: MocQty<T>> Clone for MocCell<T, Q> {
  fn clone(&self) -> MocCell<T, Q> {
    MocCell(self.0.clone(),  PhantomData)
  }
}*/


