use num::One;

use std::fmt::Debug;
use std::ops::Range;
use std::marker::PhantomData;
use crate::ranges::Idx;

/// Number of bits reserved to code the quantity type
const N_RESERVED_BITS: u8 = 2;

/// Mask selecting all but the LSB (applying the mask <=> a = if (a is even) { a } else { a - 1 }
const TO_EVEN_MASK: u32 = !0x1;

pub trait Bounded<T> {
    fn upper_bound_exclusive() -> T;
}
impl<T, Q> Bounded<T> for Q where T: Idx, Q: MocQty<T> {
    /// The largest possible value (exclusive) for a value of type T of the quantity Q.
    fn upper_bound_exclusive() -> T {
        Self::n_cells_max()
    }
}

/// Generic constants defining a quantity that can be put in a MOC,
/// independently of it the precise integer type used to represent it.
pub trait MocableQty: Send + Sync {
    /// A simple char prefix to identify the quantity (e.g. in ASCII serialisation)
    const PREFIX: char;
    /// Dimension of the qty, i.e. number of bits needed to code a sub-cell relative index
    const DIM: u8;
    /// Number of base cells, i.e. number of cell at depth 0
    /// (usually 2^dim, but 12 in the HEALPix case)
    const N_D0_CELLS: u8;
    /// Number of bits needed to code the base cell index
    const N_D0_BITS: u8 = n_bits_to_code_from_0_to_n_exclusive(Self::N_D0_CELLS);
    /// Mask to select the bit(s) of a level > 0:
    /// * dim 1: 001
    /// * dim 2: 011
    /// * dim 3: 111
    const LEVEL_MASK: u8 = (1 << Self::DIM) - 1;

    fn shift(delta_depth: u8) -> u8 {
        Self::DIM * delta_depth
    }
}

/// Returns the number of bits needed to code `n` values, with indices
/// from 0 (inclusive) to n (exclusive).
const fn n_bits_to_code_from_0_to_n_exclusive(n: u8) -> u8 {
    let n_bits_in_u8 = u8::N_BITS as u32; // = 8
    let index_max = n - 1;
    (n_bits_in_u8 - index_max.leading_zeros()) as u8
}

/// A quantity with its exact integer representation. (Replace Bounded?)
pub trait MocQty<T>: MocableQty where T: Idx
{
    const MAX_DEPTH: u8 = (T::N_BITS - (N_RESERVED_BITS + Self::N_D0_BITS)) / Self::DIM;
    const MAX_SHIFT: u32 = (Self::DIM * Self::MAX_DEPTH) as u32;
    // const MAX_VALUE : T = Self::N_D0_CELLS).into().unsigned_shl((Self::DIM * Self::MAX_DEPTH) as u32);

    // I rename max_value in n_cells_max, I could have rename in max_value_exclusive
    // (the inlcusive max_value is the value returned by this method minus one).
    fn n_cells_max() -> T {
        let nd0: T = Self::N_D0_CELLS.into();
        nd0.unsigned_shl(Self::MAX_SHIFT)
    }

    /// Upper bound on the maximum number of depths that can be coded using `n_bits`of a MOC index.
    /// I.e., maximum possible hierarchy depth on a
    /// `len = [0, 2^(delta_depth)^dim]` => `(log(len) / log(2)) / dim = delta_depth`
    fn delta_depth_max_from_n_bits(n_bits: u8) -> u8 {
        Self::delta_depth_max_from_n_bits_unchecked(n_bits).min(Self::MAX_DEPTH)
    }

    /// Same as `delta_depth_max_from_n_bits` without cecking that the result is smaller than
    /// depth_max.
    fn delta_depth_max_from_n_bits_unchecked(n_bits: u8) -> u8 {
        n_bits >> (Self::DIM - 1)
    }

    fn delta_with_depth_max(depth: u8) -> u8 {
        Self::MAX_DEPTH - depth
    }

    fn shift_from_depth_max(depth: u8) -> u8 {
        Self::shift(Self::delta_with_depth_max(depth))
    }

    // Method from former Bounded

    #[inline(always)]
    fn get_msb(x: T) -> u32 {
        T::N_BITS as u32 - x.leading_zeros() - 1
    }

    #[inline(always)]
    fn get_lsb(x: T) -> u32 {
        x.trailing_zeros() as u32
    }

    #[inline(always)]
    fn compute_min_depth(x: T) -> u8 {
        let dd = (x.trailing_zeros() as u8 / Self::DIM).min(Self::MAX_DEPTH);
        Self::MAX_DEPTH - dd
    }

/*
    #[inline(always)]
    fn get_depth(x: T) -> u32 {
        let msb = Self::get_msb(x) & TO_EVEN_MASK;
        let depth = (msb >> 1) - 1;

        depth
    }

    fn to_uniq_hpx(depth: u8, pix: T) -> T {
        let t = T::one().unsigned_shl(2).unsigned_shl((depth as u32) << 1);
        t + pix
    }

    /// Use a sentinel bit on the lowest MSB not used to code the cell.
    fn to_uniq_using_sentinel(depth: u8, cell_index: T) -> T {
        let msb: T = T::one().unsigned_shl((depth as u32) + 1);
        msb | cell_index
    }

    #[inline(always)]
    fn uniq_with_sentinel_to_range_hpx(u: T) -> Range<T> { // uniq_to_range
        let (depth, pix) = Self::hpx_cell_depth(u);
        let tdd = ((Self::HPX_MAXDEPTH - depth) << 1) as u32;
        // The length of a range computed from a pix
        // at Self::HPX_MAXDEPTH equals to 1
        Range {
            start: pix.unsigned_shl(tdd),
            end: (pix + One::one()).unsigned_shl(tdd),
        }
    }*/

}

/// HEALPix index (either Ring or Nested)
#[derive(Clone, Debug)]
pub struct Hpx<T: Idx> (std::marker::PhantomData<T>);

impl<T: Idx> MocableQty for Hpx<T> {
    const PREFIX: char = 's';
    const DIM: u8 = 2;
    const N_D0_CELLS: u8 = 12;
}

impl<T> MocQty<T> for Hpx<T> where T: Idx { }

impl<T: Idx> Hpx<T> {
    #[inline(always)]
    pub fn from_uniq_hpx(uniq: T) -> (u8, T) { // pix_depth
        let msb = Self::get_msb(uniq) & TO_EVEN_MASK;

        let depth = (msb >> 1) - 1;
        let pix = uniq - T::one().unsigned_shl(msb);

        (depth as u8, pix)
    }

    pub fn uniq_hpx_to_range(uniq: T) -> Range<T> { // uniq_to_range
        let (depth, pix) = Self::from_uniq_hpx(uniq);
        let tdd = ((Self::MAX_DEPTH - depth) << 1) as u32;
        // The length of a range computed from a pix
        // at Self::HPX_MAXDEPTH equals to 1
        Range {
            start: pix.unsigned_shl(tdd),
            end: (pix + One::one()).unsigned_shl(tdd),
        }
    }
}


/// Time index (microsec since JD=0)
#[derive(Clone, Debug)]
pub struct Time<T: Idx> (PhantomData<T>);
impl<T: Idx> MocableQty for Time<T> {
    const PREFIX: char = 't';
    const DIM: u8 = 1;
    const N_D0_CELLS: u8 = 2;
}
impl<T> MocQty<T> for Time<T> where T: Idx { }



/*pub trait Bounded<T, Q>
where
    T: MocIdx,
    Q: MocQty<T>,
{
    #[inline(always)]
    fn get_msb(x: T) -> u32 {
        num_bits::<T>() as u32 - x.leading_zeros() - 1
    }

    #[inline(always)]
    fn get_lsb(x: T) -> u32 {
        x.trailing_zeros() as u32
    }

    #[inline(always)]
    fn hpx_cell_depth(u: T) -> (u8, T) { // pix_depth
        let msb = Self::get_msb(u) & TO_EVEN_MASK;

        let depth = (msb >> 1) - 1;
        let pix = u - T::one().unsigned_shl(msb);

        (depth as u8, pix)
    }

    #[inline(always)]
    fn get_hpx_depth(x: T) -> u32 {
        let msb = Self::get_msb(x) & TO_EVEN_MASK;
        let depth = (msb >> 1) - 1;

        depth
    }

    fn to_uniq_hpx(depth: u8, pix: T) -> T {
        let t = T::one().unsigned_shl(2).unsigned_shl((depth as u32) << 1);
        t + pix
    }

    /// Use a sentinel bit on the lowest MSB not used to code the cell.
    fn to_uniq_time(depth: u8, cell_index: T) -> T {
        let msb: T = T::one().unsigned_shl((depth as u32) + 1);
        msb | cell_index
    }

    #[inline(always)]
    fn uniq_to_range_hpx(u: T) -> Range<T> { // uniq_to_range
        let (depth, pix) = Self::hpx_cell_depth(u);
        let tdd = ((Self::HPX_MAXDEPTH - depth) << 1) as u32;
        // The length of a range computed from a pix
        // at Self::HPX_MAXDEPTH equals to 1
        Range {
            start: pix.unsigned_shl(tdd),
            end: (pix + One::one()).unsigned_shl(tdd),
        }
    }
}*/


#[cfg(test)]
mod tests {
    use crate::mocqty::{MocableQty, MocQty, Time, Hpx};

    #[test]
    fn test_hpx() {
        // Independent of T
        assert_eq!(Hpx::<u64>::DIM, 2);
        assert_eq!(Hpx::<u64>::N_D0_CELLS, 12);
        assert_eq!(Hpx::<u64>::N_D0_BITS, 4);
        assert_eq!(Hpx::<u64>::LEVEL_MASK, 3);
        assert_eq!(Hpx::<u64>::shift(1), 2);
        assert_eq!(Hpx::<u64>::shift(10), 20);
        // Depends on T
        assert_eq!(Hpx::<u64>::MAX_DEPTH, 29);
        assert_eq!(Hpx::<u64>::MAX_SHIFT, 58);
        assert_eq!(Hpx::<u64>::n_cells_max(), 12 * 4_u64.pow(29));
    }

    #[test]
    fn test_time() {
        // Independent of T
        assert_eq!(Time::<u64>::DIM, 1);
        assert_eq!(Time::<u64>::N_D0_CELLS, 2);
        assert_eq!(Time::<u64>::N_D0_BITS, 1);
        assert_eq!(Time::<u64>::LEVEL_MASK, 1);
        assert_eq!(Time::<u64>::shift(1), 1);
        assert_eq!(Time::<u64>::shift(10), 10);
        // Depends on T
        assert_eq!(Time::<u64>::MAX_DEPTH, 61);
        assert_eq!(Time::<u64>::MAX_SHIFT, 61);
        assert_eq!(Time::<u64>::n_cells_max(), 2_u64.pow(62));
    }
}
