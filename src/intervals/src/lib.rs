#![feature(test)]
extern crate num;
extern crate rand;
extern crate test;
extern crate rayon;
extern crate ndarray;
extern crate healpix;

mod ranges;
mod utils;
pub mod ranges2d;
pub mod bounded;

use num::{Integer, PrimInt};
use crate::bounded::Bounded;

use std::ops::Range;
use std::mem;
use std::slice::Iter;

use rayon::prelude::*;

#[derive(Debug)]
pub struct UniqRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + Clone
    + std::fmt::Debug {
    ranges: ranges::Ranges<T>
}

impl<T> UniqRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub fn new(data: Vec<Range<T>>, min_depth: Option<i8>, make_consistent: bool) -> Self {
        let ranges = ranges::Ranges::<T>::new(data, min_depth, make_consistent);

        UniqRanges {
            ranges
        }
    }

    pub fn to_nested(self) -> NestedRanges<T> {
        let nested_data = ranges::UniqToNestedIter::new(self.ranges)
            .collect::<Vec<_>>();
        NestedRanges::<T>::new(nested_data, None, true)
    }

    pub fn iter(&self) -> Iter<Range<T>> {
        self.ranges.iter()
    }
}

impl<T> PartialEq for UniqRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn eq(&self, other: &Self) -> bool {
        self.ranges == other.ranges
    }
}

impl<T> Eq for UniqRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {}

#[derive(Debug, Clone)]
pub struct NestedRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub ranges: ranges::Ranges<T>
}

impl<T> NestedRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub fn new(data: Vec<Range<T>>, min_depth: Option<i8>, make_consistent: bool) -> Self {
        let ranges = ranges::Ranges::<T>::new(data, min_depth, make_consistent);

        NestedRanges {
            ranges
        }
    }

    pub fn to_uniq(self) -> UniqRanges<T> {
        let uniq_data = ranges::NestedToUniqIter::new(self.ranges)
            .collect::<Vec<_>>();
        UniqRanges::<T>::new(uniq_data, None, true)
    }

    pub fn depth(&self) -> i8 {
        self.ranges.depth()
    }

    pub fn degrade(&mut self, depth: i8) {
        self.ranges.degrade(depth);
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

    pub fn iter_depth_pix(self) -> ranges::DepthPixIter<T> {
        ranges::DepthPixIter::<T>::new(self.ranges)
    }

    pub fn union(&self, other: &Self) -> Self {
        let ranges = self.ranges.union(&other.ranges);
        NestedRanges {
            ranges
        }
    }

    pub fn intersection(&self, other: &Self) -> Self {
        let ranges = self.ranges.intersection(&other.ranges);
        NestedRanges {
            ranges
        }
    }

    pub fn difference(&self, other: &Self) -> Self {
        let ranges = self.ranges.difference(&other.ranges);
        NestedRanges {
            ranges
        }
    }

    pub fn complement(&self) -> Self {
        let ranges = self.ranges.complement();
        NestedRanges {
            ranges
        }
    }
}

impl<T> PartialEq for NestedRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn eq(&self, other: &Self) -> bool {
        self.ranges == other.ranges
    }
}

impl<T> Eq for NestedRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {}

pub struct RangesPy<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub data: Array2<T>,
    pub min_depth: Option<i8>,
    pub make_consistent: bool,
}

impl<T> From<RangesPy<T>> for NestedRanges<T>
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

        NestedRanges::<T>::new(data, min_depth, make_consistent)
    }
}
impl<T> From<RangesPy<T>> for UniqRanges<T>
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

        UniqRanges::<T>::new(data, min_depth, make_consistent)
    }
}

use ndarray::{Array1, Array2};
impl From<NestedRanges<u64>> for Array2<u64> {
    fn from(nested_ranges: NestedRanges<u64>) -> Self {
        let ranges = nested_ranges.ranges.0;
        // Cast Vec<Range<u64>> to Vec<u64>
        let len = ranges.len();
        let data = utils::to_flat_vec(ranges);

        // Get a Array1 from the Vec<u64> without copying any data
        let result = Array1::from_vec(data);

        // Reshape the result to get a Array2 of shape (N x 2) where N is the number 
        // of HEALPix cell contained in the moc
        result.into_shape((len, 2))
                .unwrap()
                .to_owned()
    }
}
impl From<UniqRanges<u64>> for Array2<u64> {
    fn from(uniq_ranges: UniqRanges<u64>) -> Self {
        let ranges = uniq_ranges.ranges.0;
        // Cast Vec<Range<u64>> to Vec<u64>
        let len = ranges.len();
        let data = utils::to_flat_vec(ranges);

        // Get a Array1 from the Vec<u64> without copying any data
        let result = Array1::from_vec(data);

        // Reshape the result to get a Array2 of shape (N x 2) where N is the number 
        // of HEALPix cell contained in the moc
        result.into_shape((len, 2))
                .unwrap()
                .to_owned()
    }
}

#[cfg(test)]
mod tests {
    use crate::{UniqRanges, NestedRanges};

    use num::PrimInt;

    #[test]
    fn test_uniq_iter() {
        let simple_nested = NestedRanges::<u64>::new(vec![0..1], None, true);
        let complex_nested = NestedRanges::<u64>::new(vec![7..76], None, true);
        let empty_nested = NestedRanges::<u64>::new(vec![], None, true);

        let simple_uniq = UniqRanges::<u64>::new(vec![4*4.pow(29)..(4*4.pow(29) + 1)], None, true);
        let complex_uniq = UniqRanges::<u64>::new(
            vec![
                (1 + 4*4.pow(27))..(4 + 4*4.pow(27)),
                (2 + 4*4.pow(28))..(4 + 4*4.pow(28)),
                (16 + 4*4.pow(28))..(19 + 4*4.pow(28)),
                (7 + 4*4.pow(29))..(8 + 4*4.pow(29))
            ],
            None,
            true
        );
        let empty_uniq = UniqRanges::<u64>::new(vec![], None, true);

        assert_eq!(simple_nested.clone().to_uniq(), simple_uniq);
        assert_eq!(complex_nested.clone().to_uniq(), complex_uniq);
        assert_eq!(empty_nested.clone().to_uniq(), empty_uniq);

        assert_eq!(simple_uniq.to_nested(), simple_nested);
        assert_eq!(complex_uniq.to_nested(), complex_nested);
        assert_eq!(empty_uniq.to_nested(), empty_nested);
    }

    #[test]
    fn test_uniq_nested_conversion() {
        let input = vec![1056..1057, 1057..1058, 1083..1084, 1048539..1048540, 1048574..1048575, 1048575..1048576];
        
        let ranges = UniqRanges::<u64>::new(input.clone(), None, true);
        let expected = UniqRanges::<u64>::new(input, None, true);

        assert_eq!(ranges.to_nested().to_uniq(), expected);
    }
}