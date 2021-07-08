use num::One;

use std::fmt::Debug;
use std::ops::Range;
use std::marker::PhantomData;

use crate::idx::Idx;
use crate::deser::fits::keywords::MocDim;

/// Number of bits reserved to code the quantity type
const N_RESERVED_BITS: u8 = 2;

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
pub trait MocableQty: PartialEq + Eq + Send + Sync {
    /// A simple str to identify the quantity (e.g. in ASCII serialisation)
    const NAME: &'static str;
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

    /// FITS keyword
    const MOC_DIM: MocDim;
    /// For FITS serialization (TODO: find a better approach)
    const HAS_COOSYS: bool;
    /// For FITS serialization (TODO: find a better approach)
    const HAS_TIMESYS: bool;

    /// `v * Self::DIM`, generic so that for:
    /// * `DIM=1` this is a no operation,
    /// * `DIM=2` we can use `v  << 1`
    fn mult_by_dim<T: Idx>(v: T) -> T;
    /// `v / Self::DIM`, generic so that for:
    /// * `DIM=1` this is a no operation,
    /// * `DIM=2` we can use `v  >> 1`
    fn div_by_dim<T: Idx>(v: T) -> T;

    // dim 1: delta_depth
    // dim 2: delta_depth << 1
    // dim 3:
    #[inline(always)]
    fn shift(delta_depth: u8) -> u8 {
        Self::mult_by_dim(delta_depth)
    }

}

/// Returns the number of bits needed to code `n` values, with indices
/// from 0 (inclusive) to n (exclusive).
const fn n_bits_to_code_from_0_to_n_exclusive(n: u8) -> u8 {
    let n_bits_in_u8 = u8::N_BITS as u32; // = 8
    let index_max = n - 1;
    (n_bits_in_u8 - index_max.leading_zeros()) as u8
}

/// A quantity with its exact integer representation.
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

    fn n_cells(depth: u8) -> T {
        let nd0: T = Self::N_D0_CELLS.into();
        nd0.unsigned_shl(Self::shift(depth) as u32)
    }

    /// Upper bound on the maximum number of depths that can be coded using `n_bits`of a MOC index.
    /// I.e., maximum possible hierarchy depth on a
    /// `len = [0, 2^(delta_depth)^dim]` => `(log(len) / log(2)) / dim = delta_depth`
    fn delta_depth_max_from_n_bits(n_bits: u8) -> u8 {
        Self::delta_depth_max_from_n_bits_unchecked(n_bits).min(Self::MAX_DEPTH)
    }

    /// Same as `delta_depth_max_from_n_bits` without checking that the result is smaller than
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
        let dd = Self::div_by_dim(x.trailing_zeros() as u8).min(Self::MAX_DEPTH);
        Self::MAX_DEPTH - dd
    }

    /// From generic uniq notation (using a sentinel bit)
    #[inline(always)]
    fn from_uniq_gen(uniq: T) -> (u8, T) { // pix_depth
        // T::N_BITS - uniq.leading_zeros() = number of bits to code sentinel + D + dims
        // - 1 (sentinel) - N_D0_BITS = number of bits to code dim
        let depth = Self::div_by_dim(T::N_BITS - uniq.leading_zeros() as u8 - 1 - Self::N_D0_BITS);
        let idx = uniq & !Self::sentinel_bit(depth);
        (depth as u8, idx)
    }

    /// To generic uniq notation (using a sentinel bit)
    #[inline(always)]
    fn to_uniq_gen(depth: u8, idx: T) -> T {
        Self::sentinel_bit(depth) | idx
    }

    #[inline(always)]
    fn sentinel_bit(depth: u8) -> T {
        T::one().unsigned_shl(Self::N_D0_BITS as u32).unsigned_shl(Self::shift(depth) as u32)
    }

    #[inline(always)]
    /// Range from the genric uniq notation (using a sentinel bit)
    fn uniq_gen_to_range(uniq: T) -> Range<T> { // uniq_to_range
        let (depth, pix) = Self::from_uniq_gen(uniq);
        let tdd = ((Self::MAX_DEPTH - depth) << 1) as u32;
        // The length of a range computed from a pix
        // at Self::HPX_MAXDEPTH equals to 1
        Range {
            start: pix.unsigned_shl(tdd),
            end: (pix + One::one()).unsigned_shl(tdd),
        }
    }

    /// `zuniq` is similar to the `uniq` notation (i.e. it encodes both the `depth` and
    /// the `cell index` at this depth), but the natural ordering of the type `T` preserves the
    /// global ordering of the cells, independently of the cells depth.
    /// It is similar to the [cdshealpix](https://github.com/cds-astro/cds-healpix-rust/)
    /// [BMOC](https://github.com/cds-astro/cds-healpix-rust/blob/master/src/nested/bmoc.rs)
    /// notation (without the extra bit coding a boolean)
    /// and to [multi-order-map](https://lscsoft.docs.ligo.org/ligo.skymap/moc/index.html),
    /// but also coding the depth.
    fn to_zuniq(depth: u8, idx: T) -> T {
        let zuniq = (idx << 1) | T::one();
        zuniq.unsigned_shl(Self::shift_from_depth_max(depth) as u32)
    }


    fn from_zuniq(zuniq: T) -> (u8, T) {
        let n_trailing_zero = zuniq.trailing_zeros() as u8;
        let delta_depth = Self::div_by_dim(n_trailing_zero);
        let depth = Self::MAX_DEPTH - delta_depth;
        let idx = zuniq >> (n_trailing_zero + 1) as usize;
        (depth, idx)
    }

/*
    #[inline(always)]
    fn get_depth(x: T) -> u32 {
        let msb = Self::get_msb(x) & TO_EVEN_MASK;
        let depth = (msb >> 1) - 1;

        depth
    }
*/

}

/// HEALPix index (either Ring or Nested)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hpx<T: Idx> (std::marker::PhantomData<T>);

impl<T: Idx> MocableQty for Hpx<T> {
    const NAME: &'static str = "HPX";
    const PREFIX: char = 's';
    const DIM: u8 = 2;
    const N_D0_CELLS: u8 = 12;
    // FITS specific
    const MOC_DIM: MocDim = MocDim::Space;
    const HAS_COOSYS: bool = true;
    const HAS_TIMESYS: bool = false;
    #[inline(always)]
    fn mult_by_dim<U: Idx>(v: U) -> U {
        v << 1
    }
    #[inline(always)]
    fn div_by_dim<U: Idx>(v: U) -> U {
        v >> 1
    }
}

impl<T> MocQty<T> for Hpx<T> where T: Idx { }

impl<T: Idx> Hpx<T> {

    /// From HEALPix specific uniq notation
    #[inline(always)]
    pub fn from_uniq_hpx(uniq: T) -> (u8, T) { // pix_depth
        let depth= (Self::get_msb(uniq) - 2) >> 1;
        let idx = uniq - Self::four_shl_twice_depth(depth);
        (depth as u8, idx)
    }

    /// To HEALPix specific uniq notation
    #[inline(always)]
    pub fn uniq_hpx(depth: u8, idx: T) -> T {
        idx + Self::four_shl_twice_depth(depth as u32)
    }

    #[inline(always)]
    pub fn four_shl_twice_depth(depth: u32) -> T {
        T::one().unsigned_shl(2).unsigned_shl(depth << 1)
    }

    /// Range from the HEALPix specific uniq notation
    #[inline(always)]
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
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Time<T: Idx> (PhantomData<T>);
impl<T: Idx> MocableQty for Time<T> {
    const NAME: &'static str = "TIME";
    const PREFIX: char = 't';
    const DIM: u8 = 1;
    const N_D0_CELLS: u8 = 2;
    // FITS specific
    const MOC_DIM: MocDim = MocDim::Time;
    const HAS_COOSYS: bool = false;
    const HAS_TIMESYS: bool = true;
    #[inline(always)]
    fn mult_by_dim<U: Idx>(v: U) -> U {
        v
    }
    #[inline(always)]
    fn div_by_dim<U: Idx>(v: U) -> U {
        v
    }
}
impl<T> MocQty<T> for Time<T> where T: Idx { }


#[cfg(test)]
mod tests {
    use crate::qty::{MocableQty, MocQty, Time, Hpx};

    #[test]
    fn test_hpx_uniq() {
       for depth in 0..8 {
           for idx in 0..Hpx::<u64>::n_cells(depth) {
                assert_eq!((depth, idx), Hpx::<u64>::from_uniq_hpx(Hpx::<u64>::uniq_hpx(depth, idx)));
           }
       }

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
