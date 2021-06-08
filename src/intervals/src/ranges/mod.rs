//! Very generic ranges operations

use std::{cmp, mem};
use std::str::FromStr;
use std::fmt::{Debug, Display};
use std::collections::VecDeque;
use std::ops::{Range, Index, AddAssign};
use std::slice::Iter;
use std::ptr::slice_from_raw_parts;

use num::{Integer, One, PrimInt, Zero, ToPrimitive};
use rayon::prelude::*;
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};
use ndarray::{Array1, Array2};
use byteorder::{ByteOrder, ReadBytesExt};

pub mod ranges2d;

use crate::utils;
use std::convert::TryFrom;
use std::io::Read;

// 'static mean that Idx does not contains any reference
pub trait Idx: 'static + Integer + PrimInt + ToPrimitive + AddAssign
                       + FromStr + From<u8> + TryFrom<u64>
                       + Send + Sync + Debug + Display + Copy {
    const N_BYTES: u8 = mem::size_of::<Self>() as u8;
    const N_BITS: u8 = Self::N_BYTES << 3;
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error>;
}

impl Idx for u8 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_u8()
    }
}
impl Idx for u16 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_u16::<B>()
    }
}
impl Idx for u32 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_u32::<B>()
    }
}
impl Idx for u64 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_u64::<B>()
    }
}
impl Idx for u128 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_u128::<B>()
    }
}

impl Idx for i16 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_i16::<B>()
    }
}
impl Idx for i32 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_i32::<B>()
    }
}
impl Idx for i64 {
    fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
        reader.read_i64::<B>()
    }
}

/*impl<T> Idx for T where T: 'static + Integer + PrimInt + ToPrimitive + AddAssign
                                   + FromStr + From<u8> + TryFrom<u64>
                                   + Send + Sync + Debug + Display + Copy {}*/

/// Generic operations on a set of Sorted and Non-Overlapping ranges.
/// SNO = Sorted Non-Overlapping
pub trait SNORanges<'a, T: Idx>: Sized {

    type Iter: Iterator<Item = &'a Range<T>>;
    type ParIter: ParallelIterator<Item = &'a Range<T>>;

    fn is_empty(&self) -> bool;

    /// The iterator **MUST** return sorted and non-overlapping ranges
    fn iter(&'a self) -> Self::Iter;

    fn par_iter(&'a self) -> Self::ParIter;

    fn intersects(&self, x: &Range<T>) -> bool;
    fn par_intersects(&'a self, x: &Range<T>) -> bool {
        self.par_iter()
          .map(|r| !(x.start >= r.end || x.end <= r.start))
          .any(|a| a)
    }

    fn contains_val(&self, x: &T) -> bool;

    fn contains(&self, x: &Range<T>) -> bool;
    fn par_contains(&'a self, x: &Range<T>) -> bool {
        self.par_iter()
          .map(|r| x.start >= r.start && x.end <= r.end)
          .any(|a| a)
    }

    fn merge(&self, other: &Self, op: impl Fn(bool, bool) -> bool) -> Self;

    fn union(&self, other: &Self) -> Self {
        self.merge(other, |a, b| a || b)
    }

    fn intersection(&self, other: &Self) -> Self {
        self.merge(other, |a, b| a && b)
    }

    fn difference(&self, other: &Self) -> Self {
        self.merge(other, |a, b| a && !b)
    }

    /// Returns the complement assuming that the possible values are
    /// in `[0, upper_bound_exclusive[`.
    fn complement_with_upper_bound(&self, upper_bound_exclusive: T) -> Self;

    /// Performs a logical OR between all the bounds of all ranges,
    /// and returns the number of trailing zeros.
    fn trailing_zeros(&'a self) -> u8 {
        let all_vals_or: T = self.iter()
          .fold(Zero::zero(), |res, x| res | x.start | x.end);
        all_vals_or.trailing_zeros() as u8
    }

    /// Performs a logical OR between all the bounds of all ranges,
    /// and returns the number of trailing zeros.
    fn par_trailing_zeros(&'a self) -> u8 {
        let all_vals_or: T = self.par_iter()
          .fold(|| T::zero(), |res, x| res | x.start | x.end)
          .reduce(|| T::zero(), |a, b| a | b);
        all_vals_or.trailing_zeros() as u8
    }
}

/// Generic Range operations
#[derive(Clone, Debug)]
pub struct Ranges<T: Idx>(pub Vec<Range<T>>);

impl<T: Idx> Default for Ranges<T> {
    fn default() -> Self {
        Ranges(Default::default())
    }
}

impl<T: Idx> Ranges<T> {

    /// Assumes (without checking!) that the input vector of range is already sorted and do not
    /// contains overlapping (or consecutive) ranges
    pub fn new_unchecked(data: Vec<Range<T>>) -> Self {
        Ranges(data)
    }

    /// Assumes (without checking!) that the input vector of range is already sorted **BUT**
    /// may contains overlapping (or consecutive) ranges.
    pub fn new_from_sorted(data: Vec<Range<T>>) -> Self {
        Ranges(MergeOverlappingRangesIter::new(data.iter(), None).collect::<Vec<_>>())
    }

    /// Internally sorts the input vector and ensures there is no overlapping (or consecutive) ranges.
    pub fn new_from(mut data: Vec<Range<T>>) -> Self {
        (&mut data).par_sort_unstable_by(|left, right| left.start.cmp(&right.start));
        Self::new_from_sorted(data)
    }

    /*/// Make the `Ranges<T>` consistent
    ///
    /// # Info
    ///
    /// By construction, the data are sorted so that it is possible (see the new
    /// method definition above) to merge the overlapping ranges.
    pub fn make_consistent(mut self) -> Self {
        self.0 = MergeOverlappingRangesIter::new(self.iter(), None).collect::<Vec<_>>();
        self
    }*/
}

impl<T: Idx> PartialEq for Ranges<T> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<T: Idx> Index<usize> for Ranges<T> {
    type Output = Range<T>;

    fn index(&self, index: usize) -> &Range<T> {
        &self.0[index]
    }
}

impl<'a, T: Idx> SNORanges<'a, T> for Ranges<T> {

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

    // The MOC **MUST BE** consistent!
    fn intersects(&self, x: &Range<T>) -> bool {
        let len = self.0.len() << 1;
        let ptr = self.0.as_ptr();
        let result: &[T] = unsafe {
            &* slice_from_raw_parts(ptr as *const T, len)
        };
        match result.binary_search(&x.start) {
            Ok(i) => {
                i & 1 == 0 // even index => x.start = range.start => Ok
                || (i + 1 < len && x.end > result[i + 1]) // else (odd index) => upper bound => Check x.end > next range.start
            },
            Err(i) => {
                i & 1 == 1 // index is odd => x.start inside the range
                || x.end > result[i] // else (even index) => Check x.end > range.start
            },
        }
    }

    // The MOC **MUST BE** consistent!
    fn contains_val(&self, x: &T) -> bool {
        let len = self.0.len() << 1;
        let ptr = self.0.as_ptr();
        let result: &[T] = unsafe {
            &* slice_from_raw_parts(ptr as *const T, len)
        };
        // Performs a binary search in it
        match result.binary_search(x) {
            Ok(i) => i & 1 == 0, // index must be even (lower bound of a range)
            Err(i) => i & 1 == 1, // index must be odd (inside a range)
        }
    }

    // The MOC **MUST BE** consistent!
    fn contains(&self, x: &Range<T>) -> bool {
        let len = self.0.len() << 1;
        let ptr = self.0.as_ptr();
        // It is ok because:
        // * T is a primitive type in practice, and size_of(T) << 1 == size_of(Range<T>)
        // * the alignment of Range<T> is the same as the alignment of T (see test_alignment)
        // But: https://rust-lang.github.io/unsafe-code-guidelines/layout/structs-and-tuples.html
        let result: &[T] = unsafe {
            &* slice_from_raw_parts(ptr as *const T, len)
        };
        // Performs a binary search in it
        match result.binary_search(&x.start) {
            Ok(i) => i & 1 == 0 && x.end <= result[i | 1], // index must be even (lower bound of a range)
            Err(i) => i & 1 == 1 && x.end <= result[i], // index must be odd (inside a range)
        }
    }

    fn merge(&self, other: &Self, op: impl Fn(bool, bool) -> bool) -> Self {
        let l = &self.0;
        let ll = l.len() << 1;

        let r = &other.0;
        let rl = r.len() << 1;

        let mut i = 0;
        let mut j = 0;

        let mut result = Vec::with_capacity((ll + rl) as usize);

        while i < ll || j < rl {
            let (in_l, in_r, c) = if i == ll {
                let rv = if j & 0x1 != 0 {
                    r[j >> 1].end
                } else {
                    r[j >> 1].start
                };

                let on_rising_edge_t2 = (j & 0x1) == 0;

                let in_l = false;
                let in_r = on_rising_edge_t2;

                j += 1;

                (in_l, in_r, rv)
            } else if j == rl {
                let lv = if i & 0x1 != 0 {
                    l[i >> 1].end
                } else {
                    l[i >> 1].start
                };

                let on_rising_edge_t1 = (i & 0x1) == 0;

                let in_l = on_rising_edge_t1;
                let in_r = false;

                i += 1;

                (in_l, in_r, lv)
            } else {
                let lv = if i & 0x1 != 0 {
                    l[i >> 1].end
                } else {
                    l[i >> 1].start
                };
                let rv = if j & 0x1 != 0 {
                    r[j >> 1].end
                } else {
                    r[j >> 1].start
                };

                let on_rising_edge_t1 = (i & 0x1) == 0;
                let on_rising_edge_t2 = (j & 0x1) == 0;

                let c = cmp::min(lv, rv);

                let in_l = (on_rising_edge_t1 && c == lv) | (!on_rising_edge_t1 && c < lv);
                let in_r = (on_rising_edge_t2 && c == rv) | (!on_rising_edge_t2 && c < rv);

                if c == lv {
                    i += 1;
                }
                if c == rv {
                    j += 1;
                }

                (in_l, in_r, c)
            };

            let closed = (result.len() & 0x1) == 0;

            let add = !(closed ^ op(in_l, in_r));
            if add {
                result.push(c);
            }
        }
        Ranges(utils::unflatten(&mut result))
    }

    fn complement_with_upper_bound(&self, upper_bound_exclusive: T) -> Self {
        let ranges = &self.0;
        let mut result = Vec::<Range<T>>::with_capacity((ranges.len() + 1) as usize);

        if self.is_empty() {
            result.push(T::zero()..upper_bound_exclusive);
        } else {
            let mut s = 0;
            let mut last = if ranges[0].start == T::zero() {
                s = 1;
                ranges[0].end
            } else {
                T::zero()
            };

            result = ranges.iter()
              .skip(s)
              .map(|range| {
                  let r = last..range.start;
                  last = range.end;
                  r
              })
              .collect::<Vec<_>>();

            if last < upper_bound_exclusive {
                result.push(last..upper_bound_exclusive);
            }
        }
        Ranges(result)
    }
}

impl From<Ranges<u64>> for Array2<u64> {
    fn from(input: Ranges<u64>) -> Self {
        ranges_to_array2d(input)
    }
}
impl From<Ranges<i64>> for Array2<i64> {
    fn from(input: Ranges<i64>) -> Self {
        ranges_to_array2d(input)
    }
}


pub fn ranges_to_array2d<T: Idx>(input: Ranges<T>) -> Array2<T> {
    if input.is_empty() {
        // Warning: Empty 2D numpy arrays coming from python
        // have the shape (1, 0).
        // By consistency, we also return a (1, 0) Array2 to python
        Array2::zeros((1, 0))
    } else {
        let mut ranges = input.0;
        // Cast Vec<Range<u64>> to Vec<u64>
        let len = ranges.len();
        let data = utils::flatten(&mut ranges);

        // Get a Array1 from the Vec<u64> without copying any data
        let result: Array1<T> = data.into();

        // Reshape the result to get a Array2 of shape (N x 2) where N is the number
        // of HEALPix cell contained in the moc
        result.into_shape((len, 2)).unwrap().to_owned()
    }
}


#[derive(Debug)]
pub struct MergeOverlappingRangesIter<'a, T>
where
    T: Integer,
{
    last: Option<Range<T>>,
    ranges: Iter<'a, Range<T>>,
    split_ranges: VecDeque<Range<T>>,
    shift: Option<u32>,
}

impl<'a, T> MergeOverlappingRangesIter<'a, T>
where
    T: Integer + PrimInt,
{
    pub fn new(
        mut ranges: Iter<'a, Range<T>>,
        shift: Option<u32>,
    ) -> MergeOverlappingRangesIter<'a, T> {
        let last = ranges.next().cloned();
        let split_ranges = VecDeque::<Range<T>>::new();
        MergeOverlappingRangesIter {
            last,
            ranges,
            split_ranges,
            shift,
        }
    }

    fn split_range(&self, range: Range<T>) -> VecDeque<Range<T>> {
        let mut ranges = VecDeque::<Range<T>>::new();
        match self.shift {
            None => {
                ranges.push_back(range);
            }
            Some(shift) => {
                let mut mask: T = One::one();
                mask = mask.unsigned_shl(shift) - One::one();

                if range.end - range.start < mask {
                    ranges.push_back(range);
                } else {
                    let offset = range.start & mask;
                    let mut s = range.start;
                    if offset > Zero::zero() {
                        s = (range.start - offset) + mask + One::one();
                        ranges.push_back(range.start..s);
                    }

                    while s + mask + One::one() < range.end {
                        let next = s + mask + One::one();
                        ranges.push_back(s..next);
                        s = next;
                    }

                    ranges.push_back(s..range.end);
                }
            }
        }
        ranges
    }
}

impl<'a, T> Iterator for MergeOverlappingRangesIter<'a, T>
where
    T: Integer + PrimInt,
{
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.split_ranges.is_empty() {
            return self.split_ranges.pop_front();
        }

        while let Some(curr) = self.ranges.next() {
            let prev = self.last.as_mut().unwrap();
            if curr.start <= prev.end {
                prev.end = cmp::max(curr.end, prev.end);
            } else {
                let range = self.last.clone();
                self.last = Some(curr.clone());

                self.split_ranges = self.split_range(range.unwrap());
                return self.split_ranges.pop_front();
            }
        }

        if self.last.is_some() {
            let range = self.last.clone();
            self.last = None;

            self.split_ranges = self.split_range(range.unwrap());
            return self.split_ranges.pop_front();
        } else {
            None
        }
    }
}



#[cfg(test)]
mod tests {


    use std::ops::Range;
    use std::mem::{align_of, size_of};

    use num::PrimInt;
    use rand::Rng;
    use ndarray::Array2;

    use crate::mocqty::{MocQty, Hpx};
    use crate::ranges::{SNORanges, Ranges};
    use crate::ranges::ranges_to_array2d;

    #[test]
    fn test_alignment() {
        // Ensure alignement for unsafe code in contains
        assert_eq!(align_of::<u32>(), align_of::<Range<u32>>());
        assert_eq!(align_of::<u64>(), align_of::<Range<u64>>());
        assert_eq!(align_of::<i32>(), align_of::<Range<i32>>());
        assert_eq!(align_of::<i64>(), align_of::<Range<i64>>());

        assert_eq!(size_of::<u32>() << 1, size_of::<Range<u32>>());
        assert_eq!(size_of::<u64>() << 1, size_of::<Range<u64>>());
        assert_eq!(size_of::<i32>() << 1, size_of::<Range<i32>>());
        assert_eq!(size_of::<i64>() << 1, size_of::<Range<i64>>());
    }

    #[test]
    fn empty_ranges_to_array2d() {
        let ranges = Ranges::<u64>::new_unchecked(vec![]);

        let result = ranges_to_array2d(ranges);
        assert_eq!(result, Array2::<u64>::zeros((1, 0)));
    }

    #[test]
    fn merge_range() {
        fn assert_merge(a: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let ranges = Ranges::<u64>::new_from(a);
            let expected_ranges = Ranges::<u64>::new_unchecked(expected);

            assert_eq!(ranges, expected_ranges);
        }

        assert_merge(vec![12..17, 3..5, 5..7, 6..8], vec![3..8, 12..17]);
        assert_merge(vec![0..1, 2..5], vec![0..1, 2..5]);
        assert_merge(vec![], vec![]);
        assert_merge(vec![0..6, 7..9, 8..13], vec![0..6, 7..13]);
    }


    #[test]
    fn test_union() {
        fn assert_union(a: Vec<Range<u64>>, b: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let a = Ranges::<u64>::new_from(a);
            let b = Ranges::<u64>::new_from(b);

            let expected_ranges = Ranges::<u64>::new_unchecked(expected);
            let ranges = a.union(&b);
            assert_eq!(ranges, expected_ranges);
        }

        assert_union(
            vec![12..17, 3..5, 5..7, 6..8],
            vec![0..1, 2..5],
            vec![0..1, 2..8, 12..17],
        );
        assert_union(
            vec![12..17, 3..5, 5..7, 6..8],
            vec![12..17, 3..5, 5..7, 6..8],
            vec![3..8, 12..17],
        );
        assert_union(vec![], vec![], vec![]);
        assert_union(vec![12..17], vec![], vec![12..17]);
        assert_union(vec![], vec![12..17], vec![12..17]);
        assert_union(vec![0..1, 2..3, 4..5], vec![1..22], vec![0..22]);
        assert_union(vec![0..10], vec![15..22], vec![0..10, 15..22]);
    }

    #[test]
    fn test_intersection() {
        fn assert_intersection(a: Vec<Range<u64>>, b: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let a = Ranges::<u64>::new_from(a);
            let b = Ranges::<u64>::new_from(b);

            let expected_ranges = Ranges::<u64>::new_unchecked(expected);
            let ranges = a.intersection(&b);
            assert_eq!(ranges, expected_ranges);
        }

        assert_intersection(vec![12..17, 3..5, 5..7, 6..8], vec![0..1, 2..5], vec![3..5]);
        assert_intersection(vec![], vec![0..1, 2..5], vec![]);
        assert_intersection(vec![], vec![], vec![]);
        assert_intersection(vec![2..6], vec![0..3, 4..8], vec![2..3, 4..6]);
        assert_intersection(vec![2..6], vec![2..6, 7..8], vec![2..6]);
        assert_intersection(vec![10..11], vec![10..11], vec![10..11]);
    }

    #[test]
    fn test_difference() {
        fn assert_difference(a: Vec<Range<u64>>, b: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let a = Ranges::<u64>::new_from(a);
            let b = Ranges::<u64>::new_from(b);

            let expected_ranges = Ranges::<u64>::new_unchecked(expected);
            let ranges = a.difference(&b);
            assert_eq!(ranges, expected_ranges);
        }

        assert_difference(vec![0..20], vec![5..7], vec![0..5, 7..20]);
        assert_difference(vec![0..20], vec![0..20], vec![]);
        assert_difference(vec![0..20], vec![], vec![0..20]);
        assert_difference(vec![], vec![0..5], vec![]);
        assert_difference(vec![0..20], vec![19..22], vec![0..19]);
        assert_difference(vec![0..20], vec![25..27], vec![0..20]);
        assert_difference(
            vec![0..20],
            vec![1..2, 3..4, 5..6],
            vec![0..1, 2..3, 4..5, 6..20],
        );
    }



    #[test]
    fn test_uniq_decompose() {
        macro_rules! uniq_to_pix_depth {
            ($t:ty, $size:expr) => {
                let mut rng = rand::thread_rng();

                (0..$size).for_each(|_| {
                    let depth = rng.gen_range(Range{start: 0, end: Hpx::<$t>::MAX_DEPTH});

                    let npix = 12 * 4.pow(depth as u32);
                    let pix = rng.gen_range(Range{start: 0, end: npix});

                    let uniq = 4 * 4.pow(depth as u32) + pix;
                    assert_eq!(Hpx::<$t>::from_uniq_hpx(uniq), (depth, pix));
                });
            };
        }

        uniq_to_pix_depth!(u128, 10000);
        uniq_to_pix_depth!(u64, 10000);
        uniq_to_pix_depth!(u32, 10000);
        uniq_to_pix_depth!(u8, 10000);
    }

    /*use test::Bencher;

    #[bench]
    fn bench_uniq_to_depth_pix(b: &mut Bencher) {
        let mut rng = rand::thread_rng();
        let n = test::black_box(100000);

        let uniq: Vec<u64> = (0..n)
            .map(|_| {
                let depth = rng.gen_range(0, 30);

                let npix = 12 * 4.pow(depth);
                let pix = rng.gen_range(0, npix);

                let u = 4 * 4.pow(depth) + pix;
                u
            })
            .collect();

        b.iter(|| {
            uniq.iter()
                .fold(0, |a, b| a + (u64::pix_depth(*b).0 as u64))
        });
    }*/
}
