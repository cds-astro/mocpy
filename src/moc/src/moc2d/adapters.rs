
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{
  ZSorted, NonOverlapping,
  RangeMOCIterator, CellMOCIterator, CellOrCellRangeMOCIterator,
  adapters::{
    CellMOCIteratorFromRanges,
    CellOrCellRangeMOCIteratorFromCells,
    RangeMOCIteratorFromCells,
    RangeMOCIteratorFromCellOrCellRanges
  }
};
use crate::moc2d::{
  HasTwoMaxDepth, MOC2Properties,
  RangeMOC2ElemIt, RangeMOC2Iterator, RangeMOC2IntoIterator,
  CellMOC2ElemIt, CellMOC2Iterator, CellMOC2IntoIterator,
  CellOrCellRangeMOC2ElemIt, CellOrCellRangeMOC2Iterator, CellOrCellRangeMOC2IntoIterator
};


// Range 2 Cell decorator
// - elem
pub struct RangeMOC2ElemItToCellMOC2ElemIt<
  T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
  E: RangeMOC2ElemIt<T, Q, U, R>
>(E, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>);

impl<T, Q, U, R, E> CellMOC2ElemIt<T, Q, U, R>
for RangeMOC2ElemItToCellMOC2ElemIt<T, Q, U, R, E>
  where
    T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
    E: RangeMOC2ElemIt<T, Q, U, R>
{
  type It1 = CellMOCIteratorFromRanges<T, Q, <E as RangeMOC2ElemIt<T, Q, U, R>>::It1>;
  type It2 = CellMOCIteratorFromRanges<U, R, <E as RangeMOC2ElemIt<T, Q, U, R>>::It2>;
  fn cell_mocs_it(self) -> (Self::It1, Self::It2) {
    let (it1, it2) = self.0.range_mocs_it();
    (it1.cells(), it2.cells())
  }
}

// - moc
pub struct RangeMOC2ToCellMOC2<
  T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
  K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
  L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
> (
  L, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>, PhantomData<I>, PhantomData<J>, PhantomData<K>
);
impl<T, Q, I, U, R, J, K, L> RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {
  pub fn new(it: L) -> Self {
    RangeMOC2ToCellMOC2(it, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData)
  }
}
impl<T, Q, I, U, R, J, K, L> HasTwoMaxDepth for RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  fn depth_max_1(&self) -> u8 {
    self.0.depth_max_1()
  }
  fn depth_max_2(&self) -> u8 {
    self.0.depth_max_2()
  }
}
impl<T, Q, I, U, R, J, K, L> ZSorted for RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> NonOverlapping for RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> MOC2Properties for RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> Iterator for RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  type Item = RangeMOC2ElemItToCellMOC2ElemIt<T, Q, U, R, K>;

  fn next(&mut self) -> Option<Self::Item> {
    self.0.next().map(|e|
      RangeMOC2ElemItToCellMOC2ElemIt(e, PhantomData, PhantomData, PhantomData, PhantomData)
    )
  }
}
impl<T, Q, I, U, R, J, K, L> CellMOC2Iterator<
  T, Q, CellMOCIteratorFromRanges<T, Q, <K as RangeMOC2ElemIt<T, Q, U, R>>::It1>,
  U, R, CellMOCIteratorFromRanges<U, R, <K as RangeMOC2ElemIt<T, Q, U, R>>::It2>,
  RangeMOC2ElemItToCellMOC2ElemIt<T, Q, U, R, K>
>  for RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{}

// - impl
impl<T, Q, I, U, R, J, K, L>
CellMOC2IntoIterator<
  T, Q, CellMOCIteratorFromRanges<T, Q, <K as RangeMOC2ElemIt<T, Q, U, R>>::It1>,
  U, R, CellMOCIteratorFromRanges<U, R, <K as RangeMOC2ElemIt<T, Q, U, R>>::It2>,
  RangeMOC2ElemItToCellMOC2ElemIt<T, Q, U, R, K>
> for L
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  type IntoCellMOC2Iter = RangeMOC2ToCellMOC2<T, Q, I, U, R, J, K, L>;

  fn into_cell_moc2_iter(self) -> Self::IntoCellMOC2Iter {
    RangeMOC2ToCellMOC2::new(self)
  }
}


// Range to CellOrCellRange decorator
// - elem
pub struct RangeMOC2ElemItToCellCellRangeMOC2ElemIt<
  T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
  E: RangeMOC2ElemIt<T, Q, U, R>
>(E, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>);

impl<T, Q, U, R, E> CellOrCellRangeMOC2ElemIt<T, Q, U, R>
for RangeMOC2ElemItToCellCellRangeMOC2ElemIt<T, Q, U, R, E>
  where
    T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
    E: RangeMOC2ElemIt<T, Q, U, R>
{
  type It1 = CellOrCellRangeMOCIteratorFromCells<T, Q,
    CellMOCIteratorFromRanges<T, Q, <E as RangeMOC2ElemIt<T, Q, U, R>>::It1>
  >;
  type It2 = CellOrCellRangeMOCIteratorFromCells<U, R,
    CellMOCIteratorFromRanges<U, R, <E as RangeMOC2ElemIt<T, Q, U, R>>::It2>
  >;
  fn cellcellrange_mocs_it(self) -> (Self::It1, Self::It2) {
    let (it1, it2) = self.0.range_mocs_it();
    (it1.cells().cellranges(), it2.cells().cellranges())
  }
}

// - moc
pub struct RangeMOC2ToCellCellRangeMOC2<
  T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
  K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
  L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
> (
  L, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>, PhantomData<I>, PhantomData<J>, PhantomData<K>
);
impl<T, Q, I, U, R, J, K, L> RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {
  pub fn new(it: L) -> Self {
    RangeMOC2ToCellCellRangeMOC2(it, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData)
  }
}
impl<T, Q, I, U, R, J, K, L> HasTwoMaxDepth for RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  fn depth_max_1(&self) -> u8 {
    self.0.depth_max_1()
  }
  fn depth_max_2(&self) -> u8 {
    self.0.depth_max_2()
  }
}
impl<T, Q, I, U, R, J, K, L> ZSorted for RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> NonOverlapping for RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> MOC2Properties for RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> Iterator for RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  type Item = RangeMOC2ElemItToCellCellRangeMOC2ElemIt<T, Q, U, R, K>;

  fn next(&mut self) -> Option<Self::Item> {
    self.0.next().map(|e|
      RangeMOC2ElemItToCellCellRangeMOC2ElemIt(e, PhantomData, PhantomData, PhantomData, PhantomData)
    )
  }
}
impl<T, Q, I, U, R, J, K, L> CellOrCellRangeMOC2Iterator<
  T, Q, CellOrCellRangeMOCIteratorFromCells<T, Q, CellMOCIteratorFromRanges<T, Q, <K as RangeMOC2ElemIt<T, Q, U, R>>::It1>>,
  U, R, CellOrCellRangeMOCIteratorFromCells<U, R, CellMOCIteratorFromRanges<U, R, <K as RangeMOC2ElemIt<T, Q, U, R>>::It2>>,
    RangeMOC2ElemItToCellCellRangeMOC2ElemIt<T, Q, U, R, K>
>  for RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{}

// - impl
impl<T, Q, I, U, R, J, K, L>
CellOrCellRangeMOC2IntoIterator<
  T, Q, CellOrCellRangeMOCIteratorFromCells<T, Q, CellMOCIteratorFromRanges<T, Q, <K as RangeMOC2ElemIt<T, Q, U, R>>::It1>>,
  U, R, CellOrCellRangeMOCIteratorFromCells<U, R, CellMOCIteratorFromRanges<U, R, <K as RangeMOC2ElemIt<T, Q, U, R>>::It2>>,
  RangeMOC2ElemItToCellCellRangeMOC2ElemIt<T, Q, U, R, K>
> for L
  where
    T: Idx, Q: MocQty<T>, I: RangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  type IntoCellOrCellRangeMOC2Iter = RangeMOC2ToCellCellRangeMOC2<T, Q, I, U, R, J, K, L>;

  fn into_cellcellrange_moc2_iter(self) -> Self::IntoCellOrCellRangeMOC2Iter {
    RangeMOC2ToCellCellRangeMOC2::new(self)
  }
}

// CellOrCellRange to Range decorator
// - elem
pub struct CellCellRangeMOC2ElemItToRangeMOC2ElemIt<
  T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
  E: CellOrCellRangeMOC2ElemIt<T, Q, U, R>
>(E, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>);

impl<T, Q, U, R, E> RangeMOC2ElemIt<T, Q, U, R>
for CellCellRangeMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, E>
  where
    T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
    E: CellOrCellRangeMOC2ElemIt<T, Q, U, R>
{
  type It1 = RangeMOCIteratorFromCellOrCellRanges<T, Q, <E as CellOrCellRangeMOC2ElemIt<T, Q, U, R>>::It1>;
  type It2 = RangeMOCIteratorFromCellOrCellRanges<U, R, <E as CellOrCellRangeMOC2ElemIt<T, Q, U, R>>::It2>;

  fn range_mocs_it(self) -> (Self::It1, Self::It2) {
    let (it1, it2) = self.0.cellcellrange_mocs_it();
    (it1.ranges(), it2.ranges())
  }
}

// - moc
pub struct CellCellRangeMOC2ToRangeMOC2<
  T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
  K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
  L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K>
> (
  L, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>, PhantomData<I>, PhantomData<J>, PhantomData<K>
);

impl<T, Q, I, U, R, J, K, L> CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K> {
  pub fn new(it: L) -> Self {
    CellCellRangeMOC2ToRangeMOC2(it, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData)
  }
}
impl<T, Q, I, U, R, J, K, L> HasTwoMaxDepth for CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  fn depth_max_1(&self) -> u8 {
    self.0.depth_max_1()
  }
  fn depth_max_2(&self) -> u8 {
    self.0.depth_max_2()
  }
}
impl<T, Q, I, U, R, J, K, L> ZSorted for CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> NonOverlapping for CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> MOC2Properties for CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> Iterator for CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  type Item = CellCellRangeMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, K>;

  fn next(&mut self) -> Option<Self::Item> {
    self.0.next().map(|e|
      CellCellRangeMOC2ElemItToRangeMOC2ElemIt(e, PhantomData, PhantomData, PhantomData, PhantomData)
    )
  }
}
impl<T, Q, I, U, R, J, K, L> RangeMOC2Iterator<
  T, Q, RangeMOCIteratorFromCellOrCellRanges<T, Q, <K as CellOrCellRangeMOC2ElemIt<T, Q, U, R>>::It1>,
  U, R, RangeMOCIteratorFromCellOrCellRanges<U, R, <K as CellOrCellRangeMOC2ElemIt<T, Q, U, R>>::It2>,
  CellCellRangeMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, K>
>  for CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K>
{}


// - impl
impl<T, Q, I, U, R, J, K, L>
RangeMOC2IntoIterator<
  T, Q, RangeMOCIteratorFromCellOrCellRanges<T, Q, <K as CellOrCellRangeMOC2ElemIt<T, Q, U, R>>::It1>,
  U, R, RangeMOCIteratorFromCellOrCellRanges<U, R, <K as CellOrCellRangeMOC2ElemIt<T, Q, U, R>>::It2>,
  CellCellRangeMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, K>
> for L
  where
    T: Idx, Q: MocQty<T>, I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K>
{
  type IntoRangeMOC2Iter = CellCellRangeMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>;

  fn into_range_moc2_iter(self) -> Self::IntoRangeMOC2Iter {
    CellCellRangeMOC2ToRangeMOC2::new(self)
  }
}


// Cell to Range decorator
// - elem
pub struct CellMOC2ElemItToRangeMOC2ElemIt<
  T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
  E: CellMOC2ElemIt<T, Q, U, R>
>(E, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>);

impl<T, Q, U, R, E> RangeMOC2ElemIt<T, Q, U, R>
for CellMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, E>
  where
    T: Idx, Q: MocQty<T>, U: Idx, R: MocQty<U>,
    E: CellMOC2ElemIt<T, Q, U, R>
{
  type It1 = RangeMOCIteratorFromCells<T, Q, <E as CellMOC2ElemIt<T, Q, U, R>>::It1>;
  type It2 = RangeMOCIteratorFromCells<U, R, <E as CellMOC2ElemIt<T, Q, U, R>>::It2>;

  fn range_mocs_it(self) -> (Self::It1, Self::It2) {
    let (it1, it2) = self.0.cell_mocs_it();
    (it1.ranges(), it2.ranges())
  }
}

// - moc
pub struct CellMOC2ToRangeMOC2<
  T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
  U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
  K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
  L: CellMOC2Iterator<T, Q, I, U, R, J, K>
> (
  L, PhantomData<T>, PhantomData<Q>, PhantomData<U>, PhantomData<R>, PhantomData<I>, PhantomData<J>, PhantomData<K>
);

impl<T, Q, I, U, R, J, K, L> CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K> {
  pub fn new(it: L) -> Self {
    CellMOC2ToRangeMOC2(it, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData, PhantomData)
  }
}
impl<T, Q, I, U, R, J, K, L> HasTwoMaxDepth for CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K>
{
  fn depth_max_1(&self) -> u8 {
    self.0.depth_max_1()
  }
  fn depth_max_2(&self) -> u8 {
    self.0.depth_max_2()
  }
}
impl<T, Q, I, U, R, J, K, L> ZSorted for CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> NonOverlapping for CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> MOC2Properties for CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K> {}
impl<T, Q, I, U, R, J, K, L> Iterator for CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K>
{
  type Item = CellMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, K>;

  fn next(&mut self) -> Option<Self::Item> {
    self.0.next().map(|e|
      CellMOC2ElemItToRangeMOC2ElemIt(e, PhantomData, PhantomData, PhantomData, PhantomData)
    )
  }
}
impl<T, Q, I, U, R, J, K, L> RangeMOC2Iterator<
  T, Q, RangeMOCIteratorFromCells<T, Q, <K as CellMOC2ElemIt<T, Q, U, R>>::It1>,
  U, R, RangeMOCIteratorFromCells<U, R, <K as CellMOC2ElemIt<T, Q, U, R>>::It2>,
  CellMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, K>
>  for CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K>
{}


// - impl
impl<T, Q, I, U, R, J, K, L>
RangeMOC2IntoIterator<
  T, Q, RangeMOCIteratorFromCells<T, Q, <K as CellMOC2ElemIt<T, Q, U, R>>::It1>,
  U, R, RangeMOCIteratorFromCells<U, R, <K as CellMOC2ElemIt<T, Q, U, R>>::It2>,
  CellMOC2ElemItToRangeMOC2ElemIt<T, Q, U, R, K>
> for L
  where
    T: Idx, Q: MocQty<T>, I: CellMOCIterator<T, Qty=Q>,
    U: Idx, R: MocQty<U>, J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K>
{
  type IntoRangeMOC2Iter = CellMOC2ToRangeMOC2<T, Q, I, U, R, J, K, L>;

  fn into_range_moc2_iter(self) -> Self::IntoRangeMOC2Iter {
    CellMOC2ToRangeMOC2::new(self)
  }
}


