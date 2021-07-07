
use std::marker::PhantomData;
use std::ops::Range;

use num::One;

use crate::idx::Idx;
use crate::qty::{MocQty, Hpx, Time};
use crate::elem::{
    cell::{Cell, MocCell},
    cellrange::CellRange,
    cellcellrange::{CellOrCellRange, MocCellOrCellRange}
};

// Commodity type definitions
pub type HpxRange<T> = MocRange<T, Hpx<T>>;
pub type TimeRange<T> = MocRange<T, Time<T>>;

/// Range at the deepest possible depth.
#[repr(transparent)] // To be able to transmute Vec<Range<T>> into Vec<MocRange<T, _>>
#[derive(Debug)]
pub struct MocRange<T: Idx, Q: MocQty<T>>(pub Range<T>, PhantomData<Q>);

impl<T: Idx, Q: MocQty<T>> Clone for MocRange<T, Q> {
    fn clone(&self) -> MocRange<T, Q> {
        MocRange(self.0.clone(),  PhantomData)
    }
}

impl<T: Idx, Q: MocQty<T>> MocRange<T, Q> {

    pub fn new(start: T, end: T) -> MocRange<T, Q> {
        MocRange(Range { start, end }, PhantomData)
    }

    /// # Args:
    /// * `depth_max`: the maximum depth of the cells in the range (the MocRanges `depth_max`).
    /// * `shift_dd`: bit shift of the depth difference between `Q::<T>::DEPTH_MAX` and RangeMOC `depth_max`
    /// * `range_len_min`: the length of a range of one cell (depends on the MocRanges `depth_max`).
    /// * `mask`: a mask allowing to know if the first element of a range is a cell of `depth_max`
    /// # Remark:
    /// We can deduce `shift_dd` and `range_len_min` from `Q::<T>` and `depth_max`, but we pass them
    /// to avoid recomputing them (even if it is fast)
    pub fn next_cell_with_knowledge(&mut self, depth_max: u8, shift_dd: usize, range_len_min: T, mask: T) -> Option<MocCell<T, Q>> {
        // println!("{:?}", self.0);
        let len = self.0.end - self.0.start;
        if len < T::one() {
            None
        } else if len == range_len_min || self.0.start & mask != T::zero() {
            // A range of 1 cell at depth_max
            let c = Cell::new(depth_max, self.0.start >> shift_dd).into();
            self.0.start += range_len_min;
            Some(c)
        } else {
            // dd max from number of bits to code from 0 to len
            let dd_max_from_len = Q::delta_depth_max_from_n_bits_unchecked(T::N_BITS - 1 - len.leading_zeros() as u8);
            // starting dd max from the smallest possible depth of the range lower bound
            let dd_max_from_low = Q::delta_depth_max_from_n_bits_unchecked(self.0.start.trailing_zeros() as u8);
            let delta_depth = dd_max_from_len.min(dd_max_from_low).min(Q::MAX_DEPTH);
            let delta_depth_shift = Q::shift(delta_depth) as usize;
            let c = Cell::new(Q::MAX_DEPTH - delta_depth, self.0.start >> delta_depth_shift).into();
            self.0.start += T::one() << delta_depth_shift;
            Some(c)
        }
    }


    // TODO TO BE CONTINUED...
    /*pub fn next_cellrange_with_knowledge(&mut self, depth_max: u8, shift_dd: usize, range_len_min: T, mask: T) -> Option<MocCellOrCellRange<T, Q>> {
      let len = self.0.end - self.0.start;
      if len < T::one() {
        None
      } else if len == range_len_min {
        // A range of 1 cell at depth_max
        let c = Cell::new(depth_max, self.0.start >> shift_dd).into();
        self.0.start += range_len_min;
        Some(MocCellOrCellRange::MocCell(c))
      } else if Q::DIM > 1 && self.0.start & mask != T::zero() {

      } else {
        // dd max from number of bits to code from 0 to len
        let dd_max_from_len = Q::delta_depth_max_from_n_bits_unchecked(T::N_BITS - 1 - len.leading_zeros() as u8);
        // starting dd max from the smallest possible depth of the range lower bound
        let dd_max_from_low = Q::delta_depth_max_from_n_bits_unchecked(self.0.start.trailing_zeros() as u8);
        let delta_depth = dd_max_from_len.min(dd_max_from_low).min(Q::MAX_DEPTH);
        let delta_depth_shift = Q::shift(delta_depth) as usize;



        let c = Cell::new(Q::MAX_DEPTH - delta_depth, self.0.start >> delta_depth_shift).into();
        self.0.start += T::one() << delta_depth_shift;
        Some(c)
      }
    }*/
}


impl<T: Idx, Q: MocQty<T>> From<Range<T>> for MocRange<T, Q> {
    fn from(range: Range<T>) -> Self {
        MocRange(range, PhantomData)
    }
}
impl<T: Idx, Q: MocQty<T>> From<&Range<T>> for MocRange<T, Q> {
    fn from(range: &Range<T>) -> Self {
        MocRange(range.clone(), PhantomData)
    }
}

impl<T: Idx, Q: MocQty<T>> From<(u8, T)> for MocRange<T, Q> {
    fn from((depth, ipix): (u8, T)) -> Self {
        let tdd = Q::shift_from_depth_max(depth) as u32;
        let range = Range {
            start: ipix.unsigned_shl(tdd),
            end: (ipix + One::one()).unsigned_shl(tdd),
        };
        Self::from(range)
    }
}

impl<T: Idx, Q: MocQty<T>> From<Cell<T>> for MocRange<T, Q> {
    fn from(cell: Cell<T>) -> Self {
        Self::from((cell.depth, cell.idx))
    }
}
impl<T: Idx, Q: MocQty<T>> From<&Cell<T>> for MocRange<T, Q> {
    fn from(cell: &Cell<T>) -> Self {
        Self::from((cell.depth, cell.idx))
    }
}

impl<T: Idx, Q: MocQty<T>> From<MocCell<T, Q>> for MocRange<T, Q> {
    fn from(cell: MocCell<T, Q>) -> Self {
        Self::from(cell.0)
    }
}
impl<T: Idx, Q: MocQty<T>> From<&MocCell<T, Q>> for MocRange<T, Q> {
    fn from(cell: &MocCell<T, Q>) -> Self {
        Self::from(&cell.0)
    }
}

impl<T: Idx, Q: MocQty<T>> From<CellRange<T>> for MocRange<T, Q> {
    fn from(cellrange: CellRange<T>) -> Self {
        let tdd = Q::shift_from_depth_max(cellrange.depth) as u32;
        let range = Range {
            start: cellrange.range.start.unsigned_shl(tdd),
            end: cellrange.range.end.unsigned_shl(tdd),
        };
        Self::from(range)
    }
}
impl<T: Idx, Q: MocQty<T>> From<&CellRange<T>> for MocRange<T, Q> {
    fn from(cellrange: &CellRange<T>) -> Self {
        let tdd = Q::shift_from_depth_max(cellrange.depth) as u32;
        let range = Range {
            start: cellrange.range.start.unsigned_shl(tdd),
            end: cellrange.range.end.unsigned_shl(tdd),
        };
        Self::from(range)
    }
}


impl<T: Idx, Q: MocQty<T>> From<CellOrCellRange<T>> for MocRange<T, Q> {
    fn from(cellcellrange: CellOrCellRange<T>) -> Self {
        match cellcellrange {
            CellOrCellRange::Cell(cell) => Self::from(cell),
            CellOrCellRange::CellRange(cellrange) => Self::from(cellrange),
        }
    }
}
impl<T: Idx, Q: MocQty<T>> From<&CellOrCellRange<T>> for MocRange<T, Q> {
    fn from(cellcellrange: &CellOrCellRange<T>) -> Self {
        match cellcellrange {
            CellOrCellRange::Cell(cell) => Self::from(cell),
            CellOrCellRange::CellRange(cellrange) => Self::from(cellrange),
        }
    }
}

impl<T: Idx, Q: MocQty<T>> From<MocCellOrCellRange<T, Q>> for MocRange<T, Q> {
    fn from(moccellcellrange: MocCellOrCellRange<T, Q>) -> Self {
        match moccellcellrange {
            MocCellOrCellRange::MocCell(mocell) => Self::from(mocell.0),
            MocCellOrCellRange::MocCellRange(mocellrange) => Self::from(mocellrange.0),
        }
    }
}

impl<T: Idx, Q: MocQty<T>> Iterator for MocRange<T, Q> {
    type Item = MocCell<T, Q>;

    fn next(&mut self) -> Option<Self::Item> {
        let len = self.0.end - self.0.start;
        if len < T::one() {
            None
        } else {
            // dd max from number of bits to code from 0 to len
            let dd_max_from_len = Q::delta_depth_max_from_n_bits_unchecked(T::N_BITS - 1 - len.leading_zeros() as u8);
            // starting dd max from the smallest possible depth og the range lower bound
            let dd_max_from_low = Q::delta_depth_max_from_n_bits_unchecked(self.0.start.trailing_zeros() as u8);
            let delta_depth = dd_max_from_len.min(dd_max_from_low).min(Q::MAX_DEPTH);
            let delta_depth_shift = Q::shift(delta_depth) as usize;
            let c = Cell::new(Q::MAX_DEPTH - delta_depth, self.0.start >> delta_depth_shift).into();
            self.0.start += T::one() << delta_depth_shift;
            Some(c)
        }
    }
}
