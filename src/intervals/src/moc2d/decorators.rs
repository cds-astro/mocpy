use crate::ranges::Idx;
use crate::mocqty::MocQty;
use crate::moc2d::{
  RangeMOC2ElemIt, CellMOC2ElemIt, RangeMOC2Iterator, HasTwoMaxDepth, MOC2Properties, 
  CellMOC2Iterator, CellMOC2IntoIterator, CellOrCellRangeMOC2ElemIt, CellOrCellRangeMOC2Iterator, 
  CellOrCellRangeMOC2IntoIterator};
use crate::moc::{
  CellMOCIteratorFromRanges, RangeMOCIterator, ZSorted, NonOverlapping, 
  CellOrCellRangeMOCIteratorFromCells, CellMOCIterator
};
use std::marker::PhantomData;

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



