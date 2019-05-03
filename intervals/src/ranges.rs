
use std::mem;
use std::cmp;
use std::ops::Range;
use std::slice::Iter;

use rayon::prelude::*;

use num::{Integer, PrimInt, Zero, One};
use crate::bounded::Bounded;

#[derive(Debug)]
pub struct Ranges<T>(pub Vec<Range<T>>) where T: Integer;

impl<T> Ranges<T>
where T: Integer + PrimInt + Bounded<T> + Send {
    pub fn new(mut data: Vec<Range<T>>) -> Ranges<T> {
        (&mut data).par_sort_unstable_by(|left, right| left.start.cmp(&right.start));

        let ranges: Vec<_> = MergeOverlappingRangesIter::new(data.iter()).collect();
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
    
    pub fn merge(&mut self, mut other: Self, op: &Fn(bool, bool) -> bool) {
        fn to_range_vec<T>(input: &mut Vec<T>) -> Vec<Range<T>> {
            let mut owned_input = Vec::<T>::new();
            mem::swap(&mut owned_input, input);

            let len = owned_input.len() >> 1;
            let cap = owned_input.capacity();
            let ptr = owned_input.as_mut_ptr() as *mut Range<T>;
            
            mem::forget(owned_input);

            let result = unsafe {
                Vec::from_raw_parts(ptr, len, cap)
            };
            
            result
        }
        fn to_flat_vec<T>(input: &mut Vec<Range<T>>) -> Vec<T> {
            let mut owned_input = Vec::<Range<T>>::new();
            mem::swap(&mut owned_input, input);

            let len = owned_input.len() << 1;
            let cap = owned_input.capacity();
            let ptr = owned_input.as_mut_ptr() as *mut T;
            
            mem::forget(owned_input);
            
            let result = unsafe {
                Vec::from_raw_parts(ptr, len, cap)
            };

            result
        }

        // Transmute the vectors from Vec<Range<u64>> to Vec<u64>
        let mut left = to_flat_vec(&mut self.0).into_iter();
        let mut right = to_flat_vec(&mut other.0).into_iter();

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

        self.0 = to_range_vec(&mut result);
    }

    pub fn union(&mut self, other: Self) {        
        self.merge(other, &|a, b| a || b)
    }

    pub fn intersection(&mut self, other: Self) {
        self.merge(other, &|a, b| a && b)
    }

    pub fn difference(&mut self, other: Self) {
        self.merge(other, &|a, b| a && !b)
    }

    pub fn complement(&mut self) {
        if self.is_empty() {
            self.0.push(Zero::zero()..<T>::MAXPIX);
            return;
        }

        let mut s = 0;
        let mut last = if self[0].start == Zero::zero() {
            s = 1;
            self[0].end
        } else {
            Zero::zero()
        };

        let mut result = self.0.iter().skip(s)
            .map(|range| {
                let r = last..range.start;
                last = range.end;
                r
            })
            .collect::<Vec<_>>();

        if last < <T>::MAXPIX {
            result.push(last..<T>::MAXPIX);
        }
        self.0 = result;
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> Iter<Range<T>> {
        self.0.iter()
    }
}

use std::ops::Index;
impl<T> Index<usize> for Ranges<T>
where T: Integer + PrimInt + Bounded<T> {
    type Output = Range<T>;

    fn index(&self, index: usize) -> &Range<T> {
        &self.0[index]
    }
}

use std::iter::FromIterator;
impl<T> FromIterator<Range<T>> for Ranges<T>
where T: Integer {
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
}

impl<'a, T> MergeOverlappingRangesIter<'a, T> 
where T: Integer + Clone + Copy {
    fn new(mut ranges: Iter<'a, Range<T>>) -> MergeOverlappingRangesIter<'a, T> {
        let last = ranges.next().cloned();
        MergeOverlappingRangesIter {
            last,
            ranges,
        }
    }
}

impl<'a, T> Iterator for MergeOverlappingRangesIter<'a, T> 
where T: Integer + Clone + Copy {
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(curr) = self.ranges.next() {
            let prev = self.last.as_mut().unwrap();
            if curr.start <= prev.end {
                prev.end = curr.end;
            } else {
                let next = self.last.clone();
                self.last = Some(curr.clone());
                return next
            }
        }
        if self.last.is_some() {
            let next = self.last.clone();
            self.last = None;
            next
        } else {
            None
        }
    }
}