//! Code specific to HEALPix ranges
//! Assuming the NSIDE is a power of 2, the code for NESTED indices or RING indices is the same.

use std::ops::Range;

use num::{CheckedAdd, One, Zero};

use crate::idx::Idx;
use crate::qty::{MocQty, Hpx};
use crate::elemset::range::{
    SNORanges, HpxRanges
};


pub struct HpxToUniqIter<T>
    where
        T: Idx + CheckedAdd,
{
    ranges: HpxRanges<T>,
    id: usize,
    buffer: Vec<Range<T>>,
    depth: u8,
    shift: u32,
    off: T,
    depth_off: T,
}

impl<T> HpxToUniqIter<T>
    where
        T: Idx + CheckedAdd,
{
    pub fn new(ranges: HpxRanges<T>) -> HpxToUniqIter<T> {
        let id = 0;
        let buffer = Vec::<Range<T>>::new();
        let depth = 0;
        // ((T::MAXDEPTH - depth) << 1) as u32;
        let shift = Hpx::<T>::MAX_SHIFT; // ok since depth init == 0

        let mut off: T = One::one();
        off = off.unsigned_shl(shift) - One::one();

        let mut depth_off: T = One::one();
        depth_off = depth_off.unsigned_shl(((depth + 1) << 1) as u32);

        HpxToUniqIter {
            ranges,
            id,
            buffer,

            depth,
            shift,
            off,
            depth_off,
        }
    }
}

impl<T> Iterator for HpxToUniqIter<T>
    where
        T: Idx + CheckedAdd,
{
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.ranges.is_empty() {
            let start_id = self.id;
            let end_id = self.ranges.0.0.len();
            for i in start_id..end_id {
                let range = &self.ranges.0.0[i];
                self.id += 1;

                let t1 = range.start + self.off;
                let t2 = range.end;

                let pix1 = t1.unsigned_shr(self.shift);
                let pix2 = t2.unsigned_shr(self.shift);

                let c1 = pix1.unsigned_shl(self.shift);
                let c2 = pix2.unsigned_shl(self.shift);

                if c2 > c1 {
                    self.buffer.push(c1..c2);

                    let e1 = self.depth_off.checked_add(&pix1).unwrap();
                    let e2 = self.depth_off.checked_add(&pix2).unwrap();

                    return Some(e1..e2);
                }
            }

            // new_from_sorted?
            let buffer_ranges = HpxRanges::<T>::new_from(self.buffer.clone());
            self.ranges = self.ranges.difference(&buffer_ranges);
            self.id = 0;
            self.buffer.clear();

            self.depth += 1;
            assert!(
                self.depth <= Hpx::<T>::MAX_DEPTH
                    || (self.depth > Hpx::<T>::MAX_DEPTH && self.ranges.is_empty())
            );
            if self.depth > Hpx::<T>::MAX_DEPTH && self.ranges.is_empty() {
                break;
            }

            // Recompute the constants for the new depth
            self.shift = Hpx::<T>::shift_from_depth_max(self.depth) as u32;
            self.off = One::one();
            self.off = self.off.unsigned_shl(self.shift) - One::one();

            self.depth_off = One::one();
            self.depth_off = self.depth_off.unsigned_shl(((self.depth + 1) << 1) as u32);
        }
        None
    }
}

// Iterator responsible for converting
// ranges of uniq numbers to ranges of
// nested numbers
pub struct UniqToHpxIter<T>
    where
        T: Idx + CheckedAdd,
{
    ranges: HpxRanges<T>,
    cur: T,
    id: usize,
}

impl<T> UniqToHpxIter<T>
    where
        T: Idx + CheckedAdd,
{
    pub fn new(ranges: HpxRanges<T>) -> UniqToHpxIter<T> {
        let id = 0;

        let cur = if ranges.is_empty() {
            Zero::zero()
        } else {
            ranges[id].start
        };
        UniqToHpxIter { ranges, cur, id }
    }
}

impl<T> Iterator for UniqToHpxIter<T>
    where
        T: Idx + CheckedAdd,
{
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        // Iteration through the ranges
        if self.id < self.ranges.0.0.len() {
            // We get the depth/ipix values of the
            // current uniq number
            let (depth, ipix) = Hpx::<T>::from_uniq_hpx(self.cur);

            // We compute the number of bit to shift
            let shift = ((Hpx::<T>::MAX_DEPTH - depth) << 1) as u32;

            let one: T = One::one();
            // Compute the final nested range
            // for the depth given
            let e1 = ipix.unsigned_shl(shift);
            let e2 = ipix.checked_add(&one).unwrap().unsigned_shl(shift);

            self.cur = self.cur.checked_add(&one).unwrap();

            let end = self.ranges[self.id].end;
            if self.cur == end {
                self.id += 1;

                if self.id < self.ranges.0.0.len() {
                    self.cur = self.ranges[self.id].start;
                }
            }

            Some(e1..e2)
        } else {
            None
        }
    }
}


/// IMPORTANT: the iterator is ordered according first to the cell depth
/// and then to the cell index.
/// See `ranges2cells` bench, using `CellMOCIteratorFromRanges` and then sorting the result
/// may be more efficient (x3 on the only bench done so far).
pub struct HpxUniq2DepthIdxIter<T>
    where
        T: Idx + CheckedAdd,
{
    ranges: HpxRanges<T>,
    current: Option<Range<T>>,

    last: Option<T>,
    depth: u8,
    shift: u32,

    offset: T,
    depth_offset: T,
}

impl<T> HpxUniq2DepthIdxIter<T>
    where
        T: Idx + CheckedAdd,
{
    pub fn new(ranges: HpxRanges<T>) -> HpxUniq2DepthIdxIter<T> {
        let depth = 0;
        let shift =  ((Hpx::<T>::MAX_DEPTH - depth) << 1) as u32;
        // More generic (but slower): let shift = Hpx::<T>::shift_from_depth_max(depth)

        let mut offset: T = One::one();
        offset = offset.unsigned_shl(shift) - One::one();

        let mut depth_offset: T = One::one();
        depth_offset = depth_offset.unsigned_shl(((depth + 1) << 1) as u32);

        let current = None;
        let last = None;
        HpxUniq2DepthIdxIter {
            ranges,
            current,
            last,
            depth,
            shift,
            offset,
            depth_offset,
        }
    }

    fn next_item_range(&mut self) -> Option<(i8, T)> {
        if let Some(current) = self.current.clone() {
            let last = self.last.unwrap();
            if last < current.end {
                let (depth, pix) = Hpx::<T>::from_uniq_hpx(last);
                self.last = last.checked_add(&One::one());

                Some((depth as i8, pix))
            } else {
                self.current = None;
                self.last = None;
                None
            }
        } else {
            None
        }
    }
}

impl<T> Iterator for HpxUniq2DepthIdxIter<T>
    where
        T: Idx + CheckedAdd,
{
    type Item = (i8, T);

    fn next(&mut self) -> Option<Self::Item> {
        let next_depth_pix = self.next_item_range();
        if next_depth_pix.is_some() {
            next_depth_pix
        } else {
            while !self.ranges.is_empty() {
                for range in self.ranges.iter() {
                    let t1 = range.start + self.offset;
                    let t2 = range.end;

                    let pix1 = t1.unsigned_shr(self.shift);
                    let pix2 = t2.unsigned_shr(self.shift);

                    let c1 = pix1.unsigned_shl(self.shift);
                    let c2 = pix2.unsigned_shl(self.shift);

                    if c2 > c1 {
                        let range_to_remove = HpxRanges::<T>::new_unchecked(vec![c1..c2]);
                        self.ranges = self.ranges.difference(&range_to_remove);

                        let e1 = self.depth_offset.checked_add(&pix1).unwrap();
                        let e2 = self.depth_offset.checked_add(&pix2).unwrap();

                        self.last = Some(e1);
                        self.current = Some(e1..e2);

                        return self.next_item_range();
                    }
                }
                self.depth += 1;

                // Recompute the constants for the new depth
                self.shift = ((Hpx::<T>::MAX_DEPTH - self.depth) << 1) as u32;
                // <=> let shift = Hpx::<T>::shift_from_depth_max(depth)

                self.offset = One::one();
                self.offset = self.offset.unsigned_shl(self.shift) - One::one();

                self.depth_offset = One::one();
                self.depth_offset = self.depth_offset.unsigned_shl((2 * self.depth + 2) as u32);
            }
            None
        }
    }
}


/*
// I rename Nested in Hpx because it works the same for Nested or Ring indexes.
// Since Ranges is now MocRanges and we already have the type HpxRanges
// But this is specific to HEALPix, because of methods to_uniq and complement.
#[derive(Debug, Clone)]
pub struct HpxRanges<T: MocIdx> {
    ranges: ranges::HpxRanges<T, Hpx<T>>,
}

pub type HpxRanges =
// pub struct Hpx<T: MocIdx> (std::marker::PhantomData<T>);
// pub type HpxRange<T> = MocRange<T, Hpx<T>>;
// pub type TimeRange<T> = MocRange<T, Time<T>>;

impl<T: MocIdx> HpxRanges<T> {
    pub fn new(data: Vec<Range<T>>) -> Self {
        let ranges = MocRanges::<T, Hpx<T>>::new(data);

        HpxRanges { ranges }
    }

    /// Make the NestedRanges<T> consistent
    ///
    /// # Info
    ///
    /// By construction, the data are sorted so that it is possible (see the new
    /// method definition above) to merge the overlapping ranges.
    pub fn make_consistent(mut self) -> Self {
        self.ranges = self.ranges.make_consistent();
        self
    }

    /// Divide the nested ranges into ranges of length
    /// 4**(<T>::MAXDEPTH - min_depth) if they are bigger than
    /// this size.
    ///
    /// # Info
    ///
    /// This requires min_depth to be defined between [0, <T>::MAXDEPTH]
    pub fn divide(mut self, min_depth: u8) -> Self {
        self.ranges = self.ranges.divide(min_depth);
        self
    }

    pub fn to_uniq(self) -> HpxUniqRanges<T> {
        let uniq_data = HpxToUniqIter::new(self.ranges).collect::<Vec<_>>();
        HpxUniqRanges::<T>::new(uniq_data).make_consistent()
    }

    pub fn depth(&self) -> u8 {
        self.ranges.depth()
    }

    /// # Input
    /// * `dim` dimension of the quantity
    ///     + 1 for time
    ///     + 2 for equatorial coordinates
    /// * `delta_depth`: difference between the target depth and the current depth
    pub fn degrade(&mut self, depth: u8) {
        self.ranges.degrade(depth)
    }

    pub fn intersects(&self, x: &Range<T>) -> bool {
        self.ranges.intersects(x)
    }

    pub fn contains(&self, x: &Range<T>) -> bool {
        self.ranges.contains(x)
    }

    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    pub fn iter(&self) -> Iter<Range<T>> {
        self.ranges.iter()
    }

    pub fn iter_depth_pix(self) -> HpxUniq2DepthIdxIter<T> {
        HpxUniq2DepthIdxIter::<T>::new(self.ranges)
    }

    pub fn union(&self, other: &Self) -> Self {
        let ranges = self.ranges.union(&other.ranges);
        HpxRanges { ranges }
    }

    pub fn intersection(&self, other: &Self) -> Self {
        let ranges = self.ranges.intersection(&other.ranges);
        HpxRanges { ranges }
    }

    pub fn difference(&self, other: &Self) -> Self {
        let ranges = self.ranges.difference(&other.ranges);
        HpxRanges { ranges }
    }

    pub fn complement(&self) -> Self {
        let ranges = self.ranges.complement(T::HPX_MAXPIX);
        HpxRanges { ranges }
    }
}

impl<T: MocIdx> PartialEq for HpxRanges<T> {
    fn eq(&self, other: &Self) -> bool {
        self.ranges == other.ranges
    }
}

impl<T: MocIdx> Eq for HpxRanges<T> { }

impl<T: MocIdx> From<MocRanges<T, Hpx<T>>> for HpxRanges<T> {
    fn from(ranges: MocRanges<T, Hpx<T>>) -> Self {
        HpxRanges::<T> { ranges }
    }
}

impl<T: MocIdx> From<HpxRanges<T>> for MocRanges<T, Hpx<T>> {
    fn from(hpx_ranges: HpxRanges<T>) -> Self {
        hpx_ranges.ranges
    }
}

use ndarray::Array2;

impl From<HpxRanges<u64>> for Array2<u64> {
    fn from(input: HpxRanges<u64>) -> Self {
        input.ranges.into()
    }
}
impl From<HpxRanges<i64>> for Array2<i64> {
    fn from(input: HpxRanges<i64>) -> Self {
        input.ranges.into()
    }
}
*/