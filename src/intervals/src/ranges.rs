
use std::collections::VecDeque;
use std::mem;
use std::cmp;
use std::ops::Range;
use std::slice::Iter;

use rayon::prelude::*;

use num::{One, Integer, PrimInt, Zero, CheckedAdd};
use crate::bounded::Bounded;

#[derive(Debug, Clone)]
pub struct Ranges<T>(pub Vec<Range<T>>)
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug;

impl<T> Ranges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub fn new(mut data: Vec<Range<T>>, min_depth: Option<i8>, make_consistent: bool) -> Ranges<T> {
        let ranges = if make_consistent {
            (&mut data).par_sort_unstable_by(|left, right| left.start.cmp(&right.start));

            MergeOverlappingRangesIter::new(data.iter(), min_depth).collect::<Vec<_>>()
        } else {
            data
        };
        Ranges(ranges)
    }

    pub fn to_flat_vec(self) -> Vec<T> {
        let mut data = self.0;

        let len = data.len() << 1;
        let cap = data.capacity();
        let ptr = data.as_mut_ptr() as *mut T;

        mem::forget(data);

        let result = unsafe {
            Vec::from_raw_parts(ptr, len, cap)
        };

        result
    }
    
    pub fn merge(&self, other: &Self, op: impl Fn(bool, bool) -> bool) -> Self {
        // Unflatten a stack containing T-typed elements
        // to a stack of Range<T> types elements without
        // copying the data.
        fn unflatten<T>(input: &mut Vec<T>) -> Vec<Range<T>> {
            let mut owned_input = Vec::<T>::new();
            // We swap the content refered by input with a new
            // allocated vector.
            // This fix the problem when ``input`` is freed by reaching out
            // the end of the caller scope.
            std::mem::swap(&mut owned_input, input);

            let len = owned_input.len() >> 1;
            let cap = owned_input.capacity();
            let ptr = owned_input.as_mut_ptr() as *mut Range<T>;
            
            mem::forget(owned_input);

            let result = unsafe {
                Vec::from_raw_parts(ptr, len, cap)
            };
            
            result
        }

        // Flatten a stack containing Range<T> typed elements to a stack containing
        // the start followed by the end value of the set of ranges (i.e. a Vec<T>).
        // This does a copy of the data. This is necessary because we do not want to
        // modify ``self`` as well as ``other`` and we want to return the result of 
        // the union of the two Ranges2D.
        fn flatten<T>(input: &Vec<Range<T>>) -> Vec<T>
        where T: Integer + Clone + Copy {
            input.clone()
                 .into_iter()
                 // Convert Range<T> to Vec<T> containing
                 // the start and the end values of the range.
                 .map(|r| vec![r.start, r.end])
                 // We can call flatten on a iterator containing other
                 // iterators (or collections in our case).
                 .flatten()
                 // Collect to get back a newly created Vec<T> 
                 .collect()
        }

        // Flatten the Vec<Range<u64>> to Vec<u64>.
        // This operation returns new vectors
        let mut left = flatten(&self.0).into_iter();
        let mut right = flatten(&other.0).into_iter();

        let mut left_id = 0;
        let mut right_id = 0;

        let mut result: Vec<T> = vec![];

        let mut curr_left_item = left.next();
        let mut curr_right_item = right.next();
        while curr_left_item.is_some() || curr_right_item.is_some() {
            if curr_left_item.is_some() && curr_right_item.is_some() {
                let left_item = curr_left_item.unwrap();
                let right_item = curr_right_item.unwrap();
                let curr = cmp::min(left_item, right_item);

                let in_left = !((curr < left_item) ^ ((left_id & 0x1) != 0));
                let in_right = !((curr < right_item) ^ ((right_id & 0x1) != 0));
                let in_res = op(in_left, in_right);

                if in_res ^ ((result.len() & 0x1) != 0) {
                    result.push(curr);
                }
                if curr == left_item {
                    left_id += 1;
                    curr_left_item = left.next();
                }
                if curr == right_item {
                    right_id += 1;
                    curr_right_item = right.next();
                }
            } else if curr_left_item.is_none() {
                let curr = curr_right_item.unwrap();

                let in_res = op(false, true);
                if in_res {
                    result.push(curr);
                }
                
                right_id += 1;
                curr_right_item = right.next();
            } else if curr_right_item.is_none() {
                let curr = curr_left_item.unwrap();
                
                let in_res = op(true, false);
                if in_res {
                    result.push(curr);
                }
                
                left_id += 1;
                curr_left_item = left.next();
            } else {
                unreachable!();
            }
        }

        Ranges(unflatten(&mut result))
    }

    pub fn union(&self, other: &Self) -> Self {        
        self.merge(other, |a, b| a || b)
    }

    pub fn intersection(&self, other: &Self) -> Self {
        self.merge(other, |a, b| a && b)
    }

    pub fn difference(&self, other: &Self) -> Self {
        self.merge(other, |a, b| a && !b)
    }

    pub fn complement(&self) -> Self {
        let mut result = Vec::<Range<T>>::with_capacity((self.0.len() + 1) as usize);

        if self.is_empty() {
            result.push(Zero::zero()..<T>::MAXPIX);
        } else {
            let mut s = 0;
            let mut last = if self[0].start == Zero::zero() {
                s = 1;
                self[0].end
            } else {
                Zero::zero()
            };

            result = self.0.iter()
                .skip(s)
                .map(|range| {
                    let r = last..range.start;
                    last = range.end;
                    r
                })
                .collect::<Vec<_>>();

            if last < <T>::MAXPIX {
                result.push(last..<T>::MAXPIX);
            }
        }
        Ranges(result)
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> Iter<Range<T>> {
        self.0.iter()
    }

    pub fn contains(&self, x: &Range<T>) -> bool {
        let result = self.0.par_iter()
            .map(|r| {
                if x.start >= r.end || x.end <= r.start {
                    false
                } else {
                    true
                }
            })
            .reduce_with(|a, b| a && b)
            .unwrap();
        result
    }

    pub fn degrade(&mut self, depth: i8) {
        let shift = ((<T>::MAXDEPTH - depth) << 1) as u32;

        let mut offset: T = One::one();
        offset = offset.unsigned_shl(shift) - One::one();

        let mut mask: T = One::one();
        mask = mask.checked_mul(&!offset).unwrap();

        let adda: T = Zero::zero();
        let mut addb: T = One::one();
        addb = addb.checked_mul(&offset).unwrap();

        let capacity = self.0.len();
        let mut result = Vec::<Range<T>>::with_capacity(capacity);

        for range in self.iter() {
            let a: T = range.start.checked_add(&adda).unwrap() & mask;
            let b: T = range.end.checked_add(&addb).unwrap() & mask;

            if b > a {
                result.push(a..b);
            }
        }

        self.0 = result;
    }
}

pub struct RangesPy<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub data: Array2<T>,
    pub min_depth: Option<i8>,
    pub make_consistent: bool,
}

impl<T> From<RangesPy<T>> for Ranges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn from(ranges: RangesPy<T>) -> Self {
        let mut data = ranges.data;
        let min_depth = ranges.min_depth;
        let make_consistent = ranges.make_consistent;
        
        let len = data.shape()[0];
        let cap = len;
        let ptr = data.as_mut_ptr() as *mut Range<T>;

        mem::forget(data);

        let data = unsafe {
            Vec::from_raw_parts(ptr, len, cap)
        };

        Ranges::<T>::new(data, min_depth, make_consistent)
    }
}

use ndarray::{Array1, Array2};
impl From<Ranges<u64>> for Array2<u64> {
    fn from(ranges: Ranges<u64>) -> Self {
        // Cast Vec<Range<u64>> to Vec<u64>
        let len = ranges.0.len();
        let data = ranges.to_flat_vec();

        // Get a Array1 from the Vec<u64> without copying any data
        let result = Array1::from_vec(data);

        // Reshape the result to get a Array2 of shape (N x 2) where N is the number 
        // of HEALPix cell contained in the moc
        result.into_shape((len, 2))
                .unwrap()
                .to_owned()
    }
}

impl<T> PartialEq for Ranges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

use std::ops::Index;
impl<T> Index<usize> for Ranges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    type Output = Range<T>;

    fn index(&self, index: usize) -> &Range<T> {
        &self.0[index]
    }
}

use std::iter::FromIterator;
impl<T> FromIterator<Range<T>> for Ranges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn from_iter<I: IntoIterator<Item = Range<T>>>(iter: I) -> Self {
        let mut ranges = Ranges(Vec::<Range<T>>::new());

        for range in iter {
            ranges.0.push(range);
        }

        ranges
    }
}

#[derive(Debug)]
pub struct MergeOverlappingRangesIter<'a, T>
where T: Integer + Clone + Copy {
    last: Option<Range<T>>,
    ranges: Iter<'a, Range<T>>,
    split_ranges: VecDeque<Range<T>>,
    min_depth: Option<i8>,
}

impl<'a, T> MergeOverlappingRangesIter<'a, T> 
where T: Integer + PrimInt + Clone + Copy {
    fn new(mut ranges: Iter<'a, Range<T>>, min_depth: Option<i8>) -> MergeOverlappingRangesIter<'a, T> {
        let last = ranges.next().cloned();
        let split_ranges = VecDeque::<Range<T>>::new();
	    MergeOverlappingRangesIter {
            last,
            ranges,
	        split_ranges,
	        min_depth,
        }
    }

    fn split_range(&self, range: Range<T>) -> VecDeque<Range<T>> {
    	let mut ranges = VecDeque::<Range<T>>::new();
	    match self.min_depth {
            None => { ranges.push_back(range); },
            Some(ref val) => {
                let shift = 2 * (29 - val) as u32;

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

    /*fn merge(ranges: &mut [Range<T>], idx: usize) {
        if ranges.len() > 1 {
            let m_index = (v.len() >> 1) as usize;
            let mut (l_ranges, r_ranges) = v.split_at_mut(m_index);
            rayon::join(|| merge(l_ranges),
                        || merge(r_ranges));

            // Ranges are supposed to be sorted here
            let l_index = (l_ranges.len() - 1) as usize;
            let r_index = 0 as usize;

            if l_ranges[l_index].end > r_ranges[r_index].start {
                r_ranges[r_index].start = l_ranges[l_index].start;

                ranges.swap();
            }
        }
    }*/
}

impl<'a, T> Iterator for MergeOverlappingRangesIter<'a, T> 
where T: Integer + PrimInt + Clone + Copy {
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

pub struct NestedToUniqIter<T>
where T: Integer + PrimInt + CheckedAdd
    + Bounded<T>
    + Sync + Send
    + std::fmt::Debug {
    ranges: Ranges<T>,
    id: usize,
    buffer: Vec<Range<T>>,
    depth: i8,
    shift: u32,
    off: T,
    depth_off: T,
}

impl<T> NestedToUniqIter<T>
where T: Integer + PrimInt + CheckedAdd 
        + Bounded<T>
        + Sync + Send
        + std::fmt::Debug {
    pub fn new(ranges: Ranges<T>) -> NestedToUniqIter<T> {
        let id = 0;
        let buffer = Vec::<Range<T>>::new();
        let depth = 0;
        let shift = ((T::MAXDEPTH - depth) << 1) as u32;

        let mut off: T = One::one();
        off = off.unsigned_shl(shift) - One::one();

        let mut depth_off: T = One::one();
        depth_off = depth_off.unsigned_shl((2 * depth + 2) as u32);

        NestedToUniqIter {
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

impl<T> Iterator for NestedToUniqIter<T>
where T: Integer + PrimInt + CheckedAdd
        + Bounded<T>
        + Send + Sync
        + std::fmt::Debug {
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.ranges.is_empty() {
            let start_id = self.id;
            let end_id = self.ranges.0.len();
            for i in start_id..end_id {
                let range = &self.ranges[i];
                let t1 = range.start + self.off;
                let t2 = range.end;

                let pix1 = t1.unsigned_shr(self.shift);
                let pix2 = t2.unsigned_shr(self.shift);

                let c1 = pix1.unsigned_shl(self.shift);
                let c2 = pix2.unsigned_shl(self.shift);

                self.id += 1;

                if c2 > c1 {
                    self.buffer.push(c1..c2);

                    let e1 = self.depth_off.checked_add(&pix1).unwrap();
                    let e2 = self.depth_off.checked_add(&pix2).unwrap();

                    return Some(e1..e2);
                }
            }
            
            self.ranges = self.ranges.difference(
                &Ranges::<T>::new(
                    self.buffer.clone(),
                    None,
                    true
                )
            );
            self.id = 0;
            self.buffer.clear();

            self.depth += 1;
            assert!(self.depth <= <T>::MAXDEPTH ||
                   (self.depth > <T>::MAXDEPTH && self.ranges.is_empty()));
            if self.depth > <T>::MAXDEPTH && self.ranges.is_empty() {
                break;
            }

            // Recompute the constants for the new depth
            self.shift = ((T::MAXDEPTH - self.depth) << 1) as u32;
            self.off = One::one();
            self.off = self.off.unsigned_shl(self.shift) - One::one();

            self.depth_off = One::one();
            self.depth_off = self.depth_off.unsigned_shl((2 * self.depth + 2) as u32);
        }
        None 
    }
}

pub struct DepthPixIter<T>
where T: Integer + PrimInt + CheckedAdd
    + Bounded<T>
    + Sync + Send
    + std::fmt::Debug {
    ranges: Ranges<T>,
    current: Option<Range<T>>,

    last: Option<T>,
    depth: i8,
    shift: u32,

    offset: T,
    depth_offset: T,
}

impl<T> DepthPixIter<T>
where T: Integer + PrimInt + CheckedAdd
    + Bounded<T>
    + Sync + Send
    + std::fmt::Debug {
    pub fn new(ranges: Ranges<T>) -> DepthPixIter<T> {
        let depth = 0;
        let shift = ((T::MAXDEPTH - depth) << 1) as u32;

        let mut offset: T = One::one();
        offset = offset.unsigned_shl(shift) - One::one();

        let mut depth_offset: T = One::one();
        depth_offset = depth_offset.unsigned_shl((2 * depth + 2) as u32);

        let current = None;
        let last = None;
        DepthPixIter {
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
                let (depth, pix) = <T>::pix_depth(last);
                self.last = last
                    .checked_add(&One::one());

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

impl<T> Iterator for DepthPixIter<T>
where T: Integer + PrimInt + CheckedAdd
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
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
                        self.ranges = self.ranges.difference(
                            &Ranges::<T>::new(vec![c1..c2], None, false)
                        );

                        let e1 = self.depth_offset
                            .checked_add(&pix1)
                            .unwrap();
                        let e2 = self.depth_offset
                            .checked_add(&pix2)
                            .unwrap();
                        
                        self.last = Some(e1);
                        self.current = Some(e1..e2);

                        return self.next_item_range();
                    }
                }
                self.depth += 1;

                // Recompute the constants for the new depth
                self.shift = ((T::MAXDEPTH - self.depth) << 1) as u32;
                self.offset = One::one();
                self.offset = self.offset.unsigned_shl(self.shift) - One::one();
                
                self.depth_offset = One::one();
                self.depth_offset = self.depth_offset.unsigned_shl((2 * self.depth + 2) as u32);
            }
            None
        }
    }
}

// Iterator responsible for converting
// ranges of uniq numbers to ranges of
// nested numbers
pub struct UniqToNestedIter<T>
where T: Integer + PrimInt + CheckedAdd
    + Bounded<T>
    + Sync + Send
    + std::fmt::Debug {
    ranges: Ranges<T>,
    cur: T,
    id: usize,
}

impl<T> UniqToNestedIter<T>
where T: Integer + PrimInt + CheckedAdd
    + Bounded<T>
    + Sync + Send
    + std::fmt::Debug {
    pub fn new(ranges: Ranges<T>) -> UniqToNestedIter<T> {
        let id = 0;

        let cur = if ranges.0.len() > 0 {
            ranges[id].start
        } else {
            Zero::zero()
        };
        UniqToNestedIter {
            ranges,
            cur,
            id,
        }
    }
}

impl<T> Iterator for UniqToNestedIter<T>
where T: Integer + PrimInt + CheckedAdd
    + Bounded<T>
    + Sync + Send
    + std::fmt::Debug {
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        // Iteration through the ranges
        while self.id < self.ranges.0.len() {
            // We get the depth/ipix values of the 
            // current uniq number
            let (depth, ipix) = T::pix_depth(self.cur);

            // We compute the number of bit to shift
            let shift = (T::MAXDEPTH as u32 - depth) << 1;

            let one: T = One::one();
            // Compute the final nested range
            // for the depth given
            let e1 = ipix
                .unsigned_shl(shift);
            let e2 = ipix
                .checked_add(&one)
                .unwrap()
                .unsigned_shl(shift);

            self.cur = self.cur
                .checked_add(&one)
                .unwrap();

            let end = self.ranges[self.id].end;
            if self.cur == end {
                self.id += 1;

                if self.id < self.ranges.0.len() {
                    self.cur = self.ranges[self.id].start;
                }
            }

            return Some(e1..e2)
        }
        None
    }
}

#[derive(Debug)]
pub struct Ranges2D<T, S>
where T: Integer + Clone + Copy + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Clone + Copy + Sync + Send + std::fmt::Debug {
    // First dimension
    pub x: Vec<Range<T>>,
    // Second dimension (usually the spatial one)
    // Consecutive S values do not always refer to
    // neighbours spatial cells. Therefore there is a few chance
    // that time coverage will be the same for consecutive
    // spatial cells. So the spatial cells will not be merged
    // a lot.
    pub y: Vec<Ranges<S>>,
}

type Operation<T, S> = fn(&Ranges2D<T, S>, &Ranges2D<T, S>, bool, bool, usize, usize) -> Option<Ranges<S>>;


impl<T, S> Ranges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {

    pub fn new(mut t: Vec<Range<T>>,
               s: Vec<Ranges<S>>,
               min_depth: Option<i8>,
               make_consistent: bool) -> Result<Ranges2D<T, S>, &'static str> {
        
        if t.len() != s.len() {
            return Err("Each 1st dimension range must refer to an \
                        associated set of 2nd dim ranges. \
                        Therefore `s` and `t` must be of the same size.");
        }

        let (t_ranges, s_ranges) = if make_consistent {
            // This method is responsible for creating a 2D ranges object that have the
            // following behaviour:
            // * A range on the first dimension (containing elements of type T) refers to
            //   one and only one set of ranges of the second dimension (containing elements of type S).
            // * So to speak, ranges of the first dimension that overlap are merged only if they refer
            //   to equal set of ranges of the second dimension. Otherwise, they are kept separately.
            // * A range on the first dimension referring to an empty set of ranges is discarded.
            
            // First step, the 1st dim ranges are sorted
            (&mut t).par_sort_unstable_by(|left, right| left.start.cmp(&right.start));

            // These stacks will contain the final 1st dim set of ranges
            // and for each of them, their corresponding 2nd dim set of ranges.
            let mut res_t = Vec::<Range<T>>::with_capacity(t.len());
            let mut res_s = Vec::<Ranges<S>>::with_capacity(s.len());

            // We begin to loop over the 1st dim ranges to find the first
            // range that refer to a non null 2nd dim set of ranges.
            let mut start_idx = 0 as usize;
            // This range is stored in ``prev_t`` along with its corresponding
            // set of 2nd dim ranges ``prev_s``
            let (prev_t, prev_s): (Option<Range<T>>, Option<Ranges<S>>) = {
                let mut first_t: Option<Range<T>> = None;
                let mut first_s: Option<Ranges<S>> = None;
                for (cur_t, cur_s) in t.iter().zip(s.iter()) {
                    start_idx = start_idx + 1;
                    if !cur_s.is_empty() {
                        first_t = Some(cur_t.clone());
                        first_s = Some(cur_s.clone());
                        break;
                    }
                }
                (first_t, first_s)
            };

            // If there is at least one valid 1st dim range (i.e. one that contains
            // a non null set of 2nd dim ranges).
            if prev_s.is_some() {
                let mut prev_s = prev_s.unwrap();
                let mut prev_t = prev_t.unwrap();
                // We continue looping over the 1st dim ranges from the ``start_idx`` index.
                for (cur_t, cur_s) in t.into_iter().skip(start_idx)
                                       .zip(s.into_iter().skip(start_idx)) {
                    // We discard ranges that are not valid
                    // (i.e. those associated with an empty set of 2nd dim ranges).
                    if !cur_s.is_empty() {
                        // If ``prev_t`` and ``cur_t`` are touching to each other
                        // and the 2nd dim ranges are equal
                        if cur_t.start == prev_t.end && cur_s == prev_s {
                            prev_t.end = cur_t.end;
                        // If ``prev_t`` and ``cur_t`` are overlapping
                        } else if cur_t.start < prev_t.end {
                            // We merge the ranges only if their 2nd set of ranges
                            // are matching
                            if cur_s == prev_s {
                                prev_t.end = cmp::max(prev_t.end, cur_t.end);
                            // If they are not we can be in two cases: 
                            // 1. The second range ``cur_t`` is not contained into ``prev_t``
                            // xxxx|----
                            // --|yyyyyy
                            // Therefore we will get:
                            // xx|z|yyyy
                            // where the chain of ``z`` is the union between the first and the second
                            // set of 2nd dim ranges.
                            //
                            // 2. The second range ``cur_t`` is contained into ``prev_t``
                            // xxxxxx|---
                            // -|yy|-----
                            // In this case we need to get:
                            // x|zz|x|---
                            // where the chain of ``z`` is the union between the first and the second
                            // set of 2nd dim ranges.
                            } else {
                                let t1 = prev_t.start..cur_t.start;

                                let e1 = cmp::min(cur_t.end, prev_t.end);
                                let e2 = cmp::max(cur_t.end, prev_t.end);

                                let t2 = cur_t.start..e1;
                                let t3 = e1..e2;

                                res_t.push(t1);
                                res_t.push(t2);
                                res_t.push(t3.clone());

                                // We push ``x``
                                res_s.push(prev_s.clone());
                                // We push ``z`` aka the union between
                                // ``x`` and ``y``
                                let z = prev_s.union(&cur_s);
                                res_s.push(z);
                                
                                // Depending on whether we lie on the first
                                // or the second case, we push either ``x``
                                // or ``y``
                                if e2 == prev_t.end {
                                    // 2nd case
                                    res_s.push(prev_s.clone());
                                } else {
                                    // 1st case
                                    res_s.push(cur_s.clone());
                                    prev_s = cur_s;
                                }

                                prev_t = t3;
                            }
                        // If ``prev_t`` and ``cur_t`` are not overlapping or
                        // if they are touching but the 2nd ranges are not equal
                        } else {
                            res_t.push(prev_t);
                            res_s.push(prev_s);

                            prev_t = cur_t;
                            prev_s = cur_s;
                        }
                    }
                }

                res_t.push(prev_t);
                res_s.push(prev_s);
            }

            (res_t, res_s)
        } else {
            (t, s)
        };

        Ok(Ranges2D {
            x: t_ranges,
            y: s_ranges,
        })
    }
    
    pub fn merge(&self, other: &Self, op: Operation<T, S>) -> Ranges2D<T, S> {
        // Unflatten a stack containing T-typed elements
        // to a stack of Range<T> types elements without
        // copying the data.
        fn unflatten<T>(input: &mut Vec<T>) -> Vec<Range<T>> {
            let mut owned_input = Vec::<T>::new();
            // We swap the content refered by input with a new
            // allocated vector.
            // This fix the problem when ``input`` is freed by reaching out
            // the end of the caller scope.
            std::mem::swap(&mut owned_input, input);

            let len = owned_input.len() >> 1;
            let cap = owned_input.capacity();
            let ptr = owned_input.as_mut_ptr() as *mut Range<T>;
            
            mem::forget(owned_input);

            let result = unsafe {
                Vec::from_raw_parts(ptr, len, cap)
            };
            
            result
        }

        // Flatten a stack containing Range<T> typed elements to a stack containing
        // the start followed by the end value of the set of ranges (i.e. a Vec<T>).
        // This does a copy of the data. This is necessary because we do not want to
        // modify ``self`` as well as ``other`` and we want to return the result of 
        // the union of the two Ranges2D.
        fn flatten<T>(input: &Vec<Range<T>>) -> Vec<T>
        where T: Integer + Clone + Copy {
            input.clone()
                 .into_iter()
                 // Convert Range<T> to Vec<T> containing
                 // the start and the end values of the range.
                 .map(|r| vec![r.start, r.end])
                 // We can call flatten on a iterator containing other
                 // iterators (or collections in our case).
                 .flatten()
                 // Collect to get back a newly created Vec<T> 
                 .collect()
        }

        let sentinel = <T>::MAXPIX + One::one();
        // Get the ranges from the first dimensions of self and
        // cast them to flat vectors
        let mut t1 = flatten(&self.x);
        // Push the sentinel
        t1.push(sentinel);
        // Get the first dimension ranges from other
        let mut t2 = flatten(&other.x);
        // Push the sentinel
        t2.push(sentinel);

        let mut t_ranges: Vec<T> = Vec::<T>::new();
        let mut s_ranges: Vec<Ranges<S>> = Vec::<Ranges<S>>::new();
        let mut i = 0 as usize;
        let mut j = 0 as usize;

        // We will just need a reference to the previous
        // S Ranges because we will only compare it
        // to the current S Ranges.
        // If it is equal, then we do not need to change
        // anything. If not, we have to push the previous
        // S Ranges to the resulting S Ranges stack and set
        // its value to the current S Ranges.
        let mut prev_s: Option<&Ranges<S>> = None;

        while i < t1.len() || j < t2.len() {
            let c = cmp::min(t1[i], t2[j]);
            // If the two Ranges2D have been processed
            // then we break the loop
            if c == sentinel {
                break;
            }

            let on_rising_edge_t1 = (i & 0x1) == 0;
            let on_rising_edge_t2 = (j & 0x1) == 0;
            let in_t1 = (on_rising_edge_t1 && c == t1[i]) | (!on_rising_edge_t1 && c < t1[i]);
            let in_t2 = (on_rising_edge_t2 && c == t2[j]) | (!on_rising_edge_t2 && c < t2[j]);

            let s = op(self, other, in_t1, in_t2, i, j);

            if let Some(prev_ranges) = prev_s {
                if let Some(cur_ranges) = s {
                    if !prev_ranges.eq(&cur_ranges) {
                        t_ranges.push(c);
                        t_ranges.push(c);
                        s_ranges.push(cur_ranges);
                        prev_s = s_ranges.last();
                    }
                } else {
                    t_ranges.push(c);
                    prev_s = None;
                }
            } else {
                if let Some(cur_ranges) = s {
                    t_ranges.push(c);
                    s_ranges.push(cur_ranges);
                    prev_s = s_ranges.last();
                }
            }

            if c == t1[i] {
                i += 1;
            }
            if c == t2[j] {
                j += 1;
            }
        }

        let t_ranges = unflatten(&mut t_ranges);
        Ranges2D {
            x: t_ranges,
            y: s_ranges,
        }
    }

    fn op_union(&self, other: &Self, in_t1: bool, in_t2: bool, i: usize, j: usize) -> Option<Ranges<S>> {
        if in_t1 && in_t2 {
            let s1 = &self.y[i >> 1];
            let s2 = &other.y[j >> 1];
            Some(s1.union(s2))
        } else if !in_t1 && in_t2 {
            let s2 = &other.y[j >> 1];
            Some(s2.clone())
        } else if in_t1 && !in_t2 {
            let s1 = &self.y[i >> 1];
            Some(s1.clone())
        } else {
            None
        }
    }

    pub fn union(&self, other: &Self) -> Self {
        self.merge(other, Self::op_union)
    }

    fn op_intersection(&self, other: &Self, in_t1: bool, in_t2: bool, i: usize, j: usize) -> Option<Ranges<S>> {
        if in_t1 && in_t2 {
            let s1 = &self.y[i >> 1];
            let s2 = &other.y[j >> 1];
            Some(s1.intersection(s2))
        } else {
            None
        }
    }

    pub fn intersection(&self, other: &Self) -> Self {
        self.merge(other, Self::op_intersection)
    }

    fn op_difference(&self, other: &Self, in_t1: bool, in_t2: bool, i: usize, j: usize) -> Option<Ranges<S>> {
        if in_t1 && in_t2 {
            let s1 = &self.y[i >> 1];
            let s2 = &other.y[j >> 1];
            Some(s1.difference(s2))
        } else if in_t1 && !in_t2 {
            let s1 = &self.y[i >> 1];
            Some(s1.clone())
        } else {
            None
        }
    }

    pub fn difference(&self, other: &Self) -> Self {
        self.merge(other, Self::op_difference)
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

impl<T, S> PartialEq for Ranges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
    fn eq(&self, other: &Self) -> bool {
        // It is fast to check if the two ranges
        // have the same number of ranges towards their
        // first dimension.
        if self.x.len() != other.x.len() {
            false
        } else {
            // If they have we can check if their 1st dim ranges
            // are equal
            if self.x != other.x {
                false
            } else {
                // In the last step we must verify that
                // each 2nd dim ranges are equal
                for (s1, s2) in self.y.iter()
                           .zip(other.y.iter()) {
                    if s1 != s2 {
                        return false;
                    }
                }
                true
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bounded::Bounded;
    use crate::ranges::{Ranges, Ranges2D};

    use num::{PrimInt, Integer};
    use std::ops::Range;

    fn creating_ranges<T, S>(ranges_t: Vec<Range<T>>, ranges_s: Vec<Vec<Range<S>>>) -> Ranges2D<T, S>
    where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
          S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
        let mut vec_ranges_s = Vec::<Ranges<S>>::with_capacity(ranges_t.len());

        for range_s in ranges_s.into_iter() {
            vec_ranges_s.push(Ranges::<S>::new(range_s, None, true));
        }

        Ranges2D::new(ranges_t, vec_ranges_s, None, true).unwrap()
    }

    #[test]
    fn merge_overlapping_ranges() {
        let ranges_t = vec![0..15, 7..14, 16..17, 18..19, 19..25];
        let ranges_s = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![16..21], None, true),
            Ranges::<u64>::new(vec![16..21], None, true),
            Ranges::<u64>::new(vec![16..21, 25..26], None, true),
        ];
        let ranges_2d = Ranges2D::<u64, u64>::new(ranges_t, ranges_s, None, true).unwrap();
    }

    #[test]
    fn merge_overlapping_ranges_with_empty() {
        let ranges_t = vec![0..15, 7..14, 16..17, 18..19, 19..25];
        let ranges_s = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![], None, true),
            Ranges::<u64>::new(vec![16..21], None, true),
            Ranges::<u64>::new(vec![16..21, 25..26], None, true),
        ];
        let ranges_2d = Ranges2D::<u64, u64>::new(ranges_t, ranges_s, None, true).unwrap();
    }

    #[test]
    fn creating_empty_ranges() {
        let ranges_t = vec![0..15, 7..14];
        let ranges_s = vec![
            Ranges::<u64>::new(vec![], None, true),
            Ranges::<u64>::new(vec![], None, true),
        ];
        let ranges_2d = Ranges2D::<u64, u64>::new(ranges_t, ranges_s, None, true).unwrap();
        assert!(ranges_2d.is_empty());
    }

    #[test]
    fn union_ranges_1_3() {
        let a = creating_ranges::<u64, u64>(vec![0..10], vec![vec![16..21]]);
        let b = creating_ranges::<u64, u64>(vec![10..20], vec![vec![16..21]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(vec![0..20], vec![vec![16..21]]);
        assert_eq!(res, c);
    }
    #[test]
    fn union_ranges_1_3_bis() {
        let a = creating_ranges::<u64, u64>(vec![0..10], vec![vec![16..21]]);
        let b = creating_ranges::<u64, u64>(vec![10..20], vec![vec![16..22]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(vec![0..10, 10..20], vec![vec![16..21], vec![16..22]]);
        assert_eq!(res, c);
    }
    #[test]
    fn union_ranges_covering() {
        let a = creating_ranges::<u64, u64>(vec![0..10], vec![vec![16..21]]);
        let b = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(
            vec![0..9, 9..10, 10..20],
            vec![vec![16..21], vec![0..21], vec![0..17]]
        );
        assert_eq!(res, c);
    }
    
    #[test]
    fn empty_range_union() {
        let a = creating_ranges::<u64, u64>(vec![0..1], vec![vec![]]);
        assert!(a.is_empty());
        let b = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(
            vec![9..20],
            vec![vec![0..17]]
        );
        assert_eq!(res, c);
    }

    #[test]
    fn empty_range_union_bis() {
        let b = creating_ranges::<u64, u64>(vec![0..1], vec![vec![]]);
        assert!(b.is_empty());
        let a = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(
            vec![9..20],
            vec![vec![0..17]]
        );
        assert_eq!(res, c);
    }

    #[test]
    fn complex_union() {
        let a = creating_ranges::<u64, u64>(
            vec![0..2, 3..5, 8..9, 13..14],
            vec![vec![2..3], vec![2..3], vec![5..6], vec![7..8]]
        );
        let b = creating_ranges::<u64, u64>(
            vec![1..4, 6..7, 9..10, 11..12],
            vec![vec![0..3], vec![5..6], vec![5..6], vec![10..13]]
        );

        let result = a.union(&b);
        let expected = creating_ranges::<u64, u64>(
            vec![0..1, 1..4, 4..5, 6..7, 8..10, 11..12, 13..14],
            vec![vec![2..3], vec![0..3], vec![2..3], vec![5..6], vec![5..6], vec![10..13], vec![7..8]]
        );

        assert_eq!(expected, result);
    }
}
