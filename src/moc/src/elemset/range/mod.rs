//! Code generic to the ranges of all MOC quantities

use std::ops::{Index, Range};
use std::iter::FromIterator;
use std::slice::Iter;
use std::marker::PhantomData;

use num::{One, Zero};

use crate::idx::Idx;
use crate::qty::{Bounded, MocQty, Hpx, Time};
use crate::ranges::{Ranges, SNORanges, MergeOverlappingRangesIter};

// Commodity type definitions
pub type HpxRanges<T> = MocRanges<T, Hpx<T>>;
pub type TimeRanges<T> = MocRanges<T, Time<T>>;

pub mod hpx;
pub mod uniq;

use self::uniq::HpxUniqRanges;
use self::hpx::{HpxUniq2DepthIdxIter, HpxToUniqIter};

/// Operations specific to a given quantity
#[derive(Debug)]
pub struct MocRanges<T: Idx, Q: MocQty<T>>(pub Ranges<T>, PhantomData<Q>);

impl<T: Idx, Q: MocQty<T>> MocRanges<T, Q> {
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn ranges(&self) -> &Ranges<T> {
        &self.0
    }

    pub fn into_ranges(self) -> Ranges<T> {
        self.0
    }
}

impl<T: Idx, Q: MocQty<T>> From<Ranges<T>> for MocRanges<T, Q> {
    fn from(ranges: Ranges<T>) -> Self {
        MocRanges(ranges,  PhantomData)
    }
}

impl<T: Idx, Q: MocQty<T>> Clone for MocRanges<T, Q> {
    fn clone(&self) -> MocRanges<T, Q> {
        MocRanges(self.0.clone(),  PhantomData)
    }
}

impl<T: Idx, Q: MocQty<T>> Default for MocRanges<T, Q> {
    fn default() -> Self {
        MocRanges(Default::default(),  PhantomData)
    }
}

impl<'a, T: Idx, Q: MocQty<T>> SNORanges<'a, T> for MocRanges<T, Q> {

    type Iter = Iter<'a, Range<T>>;
    type ParIter = rayon::slice::Iter<'a, Range<T>>;

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn iter(&'a self) -> Self::Iter {
        self.0.iter()
    }

    fn par_iter(&'a self) -> Self::ParIter {
        self.0.par_iter()
    }

    fn union(&self, other: &Self) -> Self {
        self.0.union(&other.0).into()
    }

    fn intersection(&self, other: &Self) -> Self {
        self.0.intersection(&other.0).into()
    }

    fn intersects(&self, x: &Range<T>) -> bool {
        self.0.intersects(x)
    }

    fn contains_val(&self, x: &T) -> bool {
        self.0.contains_val(x)
    }

    fn contains(&self, x: &Range<T>) -> bool {
        self.0.contains(x)
    }

    fn merge(&self, other: &Self, op: impl Fn(bool, bool) -> bool) -> Self {
        Self(self.0.merge(&other.0, op), PhantomData)
    }

    fn complement_with_upper_bound(&self, upper_bound_exclusive: T) -> Self {
        Self(self.0.complement_with_upper_bound(upper_bound_exclusive), PhantomData)
    }
}


impl<T: Idx, Q: MocQty<T>> MocRanges<T, Q> {

    /// Assumes (without checking!) that the input vector of range is already sorted and do not
    /// contains overlapping (or consecutive) ranges
    pub fn new_unchecked(data: Vec<Range<T>>) -> Self {
        MocRanges(Ranges::new_unchecked(data), PhantomData)
    }

    /// Assumes (without checking!) that the input vector of range is already sorted **BUT**
    /// may contains overlapping (or consecutive) ranges.
    pub fn new_from_sorted(data: Vec<Range<T>>) -> Self {
        MocRanges(Ranges::new_from_sorted(data), PhantomData)
    }

    /// Internally sorts the input vector and ensures there is no overlapping (or consecutive) ranges.
    pub fn new_from(data: Vec<Range<T>>) -> Self {
        MocRanges(Ranges::new_from(data), PhantomData)
    }

    /// Divide the nested ranges into ranges of length
    /// `4**(<T>::MAXDEPTH - min_depth)`
    ///
    /// # Info
    ///
    /// This requires min_depth to be defined between `[0, <T>::MAXDEPTH]`
    //pub fn divide(mut self, min_depth: i8) -> Self {
    pub fn divide(mut self, min_depth: u8) -> Self {
        self.0 = Ranges::new_unchecked(
            MergeOverlappingRangesIter::new(self.iter(), Some(Q::shift_from_depth_max(min_depth) as u32))
                .collect::<Vec<_>>()
        );
        self
    }

    pub fn complement(&self) -> Self {
        self.complement_with_upper_bound(Q::upper_bound_exclusive())
    }

    pub fn iter(&self) -> Iter<Range<T>> {
        self.0.iter()
    }

    // Because of degrade, this is not a simple range but contains the notion of order/level/depth.
    pub fn degraded(&self, depth: u8) -> Self {
        let shift = Q::shift_from_depth_max(depth) as u32;

        let mut offset: T = One::one();
        offset = offset.unsigned_shl(shift) - One::one();

        let mut mask: T = One::one();
        mask = mask.checked_mul(&!offset).unwrap();

        let adda: T = Zero::zero();
        let mut addb: T = One::one();
        addb = addb.checked_mul(&offset).unwrap();

        let capacity = self.0.0.len();
        let mut result = Vec::<Range<T>>::with_capacity(capacity);

        for range in self.iter() {
            let a: T = range.start.checked_add(&adda).unwrap() & mask;
            let b: T = range.end.checked_add(&addb).unwrap() & mask;

            // if b > a {
            result.push(a..b);
            // }
        }

        // TODO: change the algo: one can merge while degrading!
        Ranges::new_unchecked(
            MergeOverlappingRangesIter::new(result.iter(), None)
              .collect::<Vec<_>>()
        ).into()
    }
    
    pub fn degrade(&mut self, depth: u8) {
        self.0 = self.degraded(depth).0
    }

    pub fn compute_min_depth(&self) -> u8 {
        Self::compute_min_depth_gen(&self.0)
        // Q::MAX_DEPTH - (self.trailing_zeros() / Q::DIM).min(Q::MAX_DEPTH)
    }

    pub fn compute_min_depth_gen(ranges: &Ranges<T>) -> u8 {
        Q::MAX_DEPTH - (ranges.trailing_zeros() / Q::DIM).min(Q::MAX_DEPTH)
    }

}

impl<T: Idx> MocRanges<T, Hpx<T>> {
    /*pub fn to_hpx_uniq_iter(self) -> HpxToUniqIter<T> {
        HpxToUniqIter::new(self.ranges).collect::<Vec<_>>();
        HpxUniqRanges::<T>::new(uniq_data).make_consistent()
    }*/
    pub fn into_hpx_uniq(self) -> HpxUniqRanges<T>{
        let uniq_data = HpxToUniqIter::new(self).collect::<Vec<_>>();
        HpxUniqRanges::<T>::new_from_sorted(uniq_data)
    }
    pub fn iter_depth_pix(self) -> HpxUniq2DepthIdxIter<T> {
        HpxUniq2DepthIdxIter::<T>::new(self)
    }
}


impl<T: Idx, Q: MocQty<T>> PartialEq for MocRanges<T, Q> {
    fn eq(&self, other: &Self) -> bool {
        self.0.0 == other.0.0
    }
}

impl<T: Idx, Q: MocQty<T>> Index<usize> for MocRanges<T, Q> {
    type Output = Range<T>;

    fn index(&self, index: usize) -> &Range<T> {
        &self.0.0[index]
    }
}

impl<T, Q> FromIterator<Range<T>> for MocRanges<T, Q>
    where
        T: Idx,
        Q: MocQty<T>{
    fn from_iter<I: IntoIterator<Item = Range<T>>>(iter: I) -> Self {
        MocRanges(
            Ranges::new_unchecked(iter.into_iter().collect::<Vec<Range<T>>>()),
            PhantomData
        )
    }
}

/*
impl<Q: MocQty<u64>> From<MocRanges<u64, Q>> for Array2<u64> {
  fn from(input: MocRanges<u64, Q>) -> Self {
    ranges_to_array2d(input.0)
  }
}
impl<Q: MocQty<i64>> From<MocRanges<i64, Q>> for Array2<i64> {
  fn from(input: MocRanges<i64, Q>) -> Self {
    ranges_to_array2d(input.0)
  }
}
*/



/*
pub struct MocRangesRef<'a, T: Idx, Q: MocQty<T>>(pub &'a Ranges<T>, PhantomData<Q>);

impl<'a, T: Idx, Q: MocQty<T>> From<&'a Ranges<T>> for MocRangesRef<'a, T, Q> {
  fn from(ranges: &'a Ranges<T>) -> Self {
    MocRangesRef(ranges,  PhantomData)
  }
}
*/



#[cfg(test)]
mod tests {
    use std::ops::Range;

    use crate::qty::{Hpx, MocQty};
    use crate::elemset::range::HpxRanges;

    #[test]
    fn merge_range_min_depth() {
        let ranges = HpxRanges::<u64>::new_unchecked(vec![0..(1 << 58)])
            .divide(1);
        let expected_ranges = vec![
            0..(1 << 56),
            (1 << 56)..(1 << 57),
            (1 << 57)..3 * (1 << 56),
            3 * (1 << 56)..(1 << 58),
        ];

        assert_eq!(ranges.0.0, expected_ranges.into_boxed_slice());
    }

    #[test]
    fn test_complement() {
        fn assert_complement(input: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let ranges = HpxRanges::<u64>::new_from_sorted(input);
            let expected_ranges = HpxRanges::<u64>::new_from_sorted(expected);

            let result = ranges.complement();
            assert_eq!(result, expected_ranges);
        }

        fn assert_complement_pow_2(input: Vec<Range<u64>>) {
            let ranges = HpxRanges::<u64>::new_from_sorted(input.clone());
            let start_ranges = HpxRanges::<u64>::new_from_sorted(input);

            let result = ranges.complement();
            let result = result.complement();

            assert_eq!(result, start_ranges);
        }

        assert_complement(vec![5..10], vec![0..5, 10..Hpx::<u64>::n_cells_max()]);
        assert_complement_pow_2(vec![5..10]);

        assert_complement(vec![0..10], vec![10..Hpx::<u64>::n_cells_max()]);
        assert_complement_pow_2(vec![0..10]);

        assert_complement(
            vec![0..1, 2..3, 4..5, 6..Hpx::<u64>::n_cells_max()],
            vec![1..2, 3..4, 5..6],
        );
        assert_complement_pow_2(vec![0..1, 2..3, 4..5, 6..Hpx::<u64>::n_cells_max()]);

        assert_complement(vec![], vec![0..Hpx::<u64>::n_cells_max()]);
        assert_complement_pow_2(vec![]);

        assert_complement(vec![0..Hpx::<u64>::n_cells_max()], vec![]);
    }

    #[test]
    fn test_depth() {
        let r1 = HpxRanges::<u64>::new_unchecked(vec![0_u64..4 * 4_u64.pow(29 - 1)]);
        assert_eq!(r1.compute_min_depth(), 0);

        let r2 = HpxRanges::<u64>::new_unchecked(vec![0_u64..4 * 4_u64.pow(29 - 3)]);
        assert_eq!(r2.compute_min_depth(), 2);

        let r3 = HpxRanges::<u64>::new_unchecked(vec![0_u64..3 * 4_u64.pow(29 - 3)]);
        assert_eq!(r3.compute_min_depth(), 3);

        let r4 = HpxRanges::<u64>::new_unchecked(vec![0_u64..12 * 4_u64.pow(29)]);
        assert_eq!(r4.compute_min_depth(), 0);

        let r5 = HpxRanges::<u64>::default();
        assert_eq!(r5.compute_min_depth(), 0);
    }

    #[test]
    fn test_degrade() {
        let mut r1 = HpxRanges::<u64>::new_unchecked(vec![0_u64..4 * 4_u64.pow(29 - 1)]);
        r1.degrade(0);
        assert_eq!(r1.compute_min_depth(), 0);

        let mut r2 = HpxRanges::<u64>::new_unchecked(vec![0_u64..4 * 4_u64.pow(29 - 3)]);
        r2.degrade(1);
        assert_eq!(r2.compute_min_depth(), 1);

        let mut r3 = HpxRanges::<u64>::new_unchecked(vec![0_u64..3 * 4_u64.pow(29 - 3)]);
        r3.degrade(1);
        assert_eq!(r3.compute_min_depth(), 1);

        let mut r4 = HpxRanges::<u64>::new_unchecked(vec![0_u64..12 * 4_u64.pow(29)]);
        r4.degrade(0);
        assert_eq!(r4.compute_min_depth(), 0);

        let mut r5 = HpxRanges::<u64>::new_unchecked(vec![0_u64..4 * 4_u64.pow(29 - 3)]);
        r5.degrade(5);
        assert_eq!(r5.compute_min_depth(), 2);
    }
}
