
use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};
use std::marker::PhantomData;

/// Decorator that converts from one datatype to another datatype.
/*pub fn convert_to_u64<T, Q, I, R>(it: I) -> ConvertIterator<T, Q, I, u64, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{
  ConvertIterator::<T, Q, I, u64, R>::new(it)
}*/

/// Decorator that converts from one datatype to another datatype.
pub fn convert<T, Q, I, U, R>(it: I) -> ConvertIterator<T, Q, I, U, R>
  where
      T: Idx,
      Q: MocQty<T>,
      I: RangeMOCIterator<T, Qty=Q>,
      U: Idx + From<T>,
      R: MocQty<U>,
{
  ConvertIterator::new(it)
}

/// Decorator that converts from one datatype to another datatype.
pub struct ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{
  last: Option<Range<U>>,
  it: I,
  _t_type: PhantomData<T>,
  _q_type: PhantomData<Q>,
  _r_type: PhantomData<R>,
}

impl<T, Q, I, U, R> ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{

  pub fn new(it: I) -> ConvertIterator<T, Q, I, U, R>  {
    let last = it.peek_last().map(|r| r.start.convert()..r.end.convert());
    ConvertIterator {
      last,
      it,
      _t_type: PhantomData,
      _q_type: PhantomData,
      _r_type: PhantomData
    }
  }
}

impl<T, Q, I, U, R> HasMaxDepth for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T, Q, I, U, R> ZSorted for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{ }
impl<T, Q, I, U, R> NonOverlapping for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{ }
impl<T, Q, I, U, R> MOCProperties for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{ }
impl<T, Q, I, U, R> Iterator for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{

  type Item = Range<U>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next().map(|r| r.start.convert()..r.end.convert())
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    self.it.size_hint()
  }
}

impl<T, Q, I, U, R> RangeMOCIterator<U> for ConvertIterator<T, Q, I, U, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx + From<T>,
    R: MocQty<U>,
{
  type Qty = R;

  fn peek_last(&self) -> Option<&Range<U>> {
    self.last.as_ref()
  }
}


/// Decorator that converts from u64 to another datatype.
pub fn convert_from_u64<Q1, T, Q, I>(it: I) -> ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{
  ConvertFromU64Iterator::new(it)
}

/// Decorator that converts from one datatype to another datatype.
pub struct ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{
  last: Option<Range<T>>,
  it: I,
  _t_type: PhantomData<T>,
  _q_type: PhantomData<Q>,
  _q1_type: PhantomData<Q1>,
}

impl<Q1, T, Q, I> ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{

  pub fn new(it: I) -> ConvertFromU64Iterator<Q1, T, Q, I>  {
    let last = it.peek_last().map(|r| T::from_u64_idx(r.start)..T::from_u64_idx(r.end));
    ConvertFromU64Iterator {
      last,
      it,
      _t_type: PhantomData,
      _q_type: PhantomData,
      _q1_type: PhantomData
    }
  }
}

impl<Q1, T, Q, I> HasMaxDepth for ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<Q1, T, Q, I> ZSorted for ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{ }
impl<Q1, T, Q, I> NonOverlapping for ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{ }
impl<Q1, T, Q, I> MOCProperties for ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{ }
impl<Q1, T, Q, I> Iterator for ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{

  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next().map(|r| T::from_u64_idx(r.start)..T::from_u64_idx(r.end))
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    self.it.size_hint()
  }
}

impl<Q1, T, Q, I> RangeMOCIterator<T> for ConvertFromU64Iterator<Q1, T, Q, I>
  where
    Q1: MocQty<u64>,
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<u64, Qty=Q1>,
{
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    self.last.as_ref()
  }
}


/// Decorator that converts from u64 to another datatype.
pub fn convert_to_u64<T, Q, I, R>(it: I) -> ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{
  ConvertToU64Iterator::new(it)
}

/// Decorator that converts from one datatype to another datatype.
pub struct  ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{
  last: Option<Range<u64>>,
  it: I,
  _t_type: PhantomData<T>,
  _q_type: PhantomData<Q>,
  _r_type: PhantomData<R>,
}

impl<T, Q, I, R> ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{

  pub fn new(it: I) -> ConvertToU64Iterator<T, Q, I, R>  {
    let last = it.peek_last().map(|r| r.start.to_u64_idx()..r.end.to_u64_idx());
    ConvertToU64Iterator {
      last,
      it,
      _t_type: PhantomData,
      _q_type: PhantomData,
      _r_type: PhantomData
    }
  }
}

impl<T, Q, I, R> HasMaxDepth for  ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}
impl<T, Q, I, R> ZSorted for  ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{ }
impl<T, Q, I, R> NonOverlapping for  ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{ }
impl<T, Q, I, R> MOCProperties for  ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{ }
impl<T, Q, I, R> Iterator for  ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{

  type Item = Range<u64>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next().map(|r| r.start.to_u64_idx()..r.end.to_u64_idx())
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    self.it.size_hint()
  }
}

impl<T, Q, I, R> RangeMOCIterator<u64> for  ConvertToU64Iterator<T, Q, I, R>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    R: MocQty<u64>,
{
  type Qty = R;

  fn peek_last(&self) -> Option<&Range<u64>> {
    self.last.as_ref()
  }
}
