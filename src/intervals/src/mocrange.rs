
use std::marker::PhantomData;
use std::ops::Range;

use num::One;

use crate::mocqty::{MocQty, Hpx, Time};
use crate::ranges::Idx;

// Commodity type definitions
pub type HpxRange<T> = MocRange<T, Hpx<T>>;
pub type TimeRange<T> = MocRange<T, Time<T>>;

#[derive(Debug)]
pub struct MocRange<T: Idx, Q: MocQty<T>>(pub Range<T>, PhantomData<Q>);
impl<T: Idx, Q: MocQty<T>> Clone for MocRange<T, Q> {
  fn clone(&self) -> MocRange<T, Q> {
    MocRange(self.0.clone(),  PhantomData)
  }
}

impl<T: Idx, Q: MocQty<T>> From<(u8, T)> for MocRange<T, Q> {
  fn from((depth, ipix): (u8, T)) -> Self {
    let tdd = Q::shift_from_depth_max(depth) as u32;
    MocRange(
      Range {
        start: ipix.unsigned_shl(tdd),
        end: (ipix + One::one()).unsigned_shl(tdd),
      },
      PhantomData
    )
  }
}
