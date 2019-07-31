use num::{Integer, PrimInt};

use crate::ranges::{Ranges, UniqToNestedIter};
use crate::nestedranges::NestedRanges;
use crate::bounded::Bounded;

use std::ops::Range;
use std::slice::Iter;

#[derive(Debug)]
pub struct UniqRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + Clone
    + std::fmt::Debug {
    ranges: Ranges<T>
}

impl<T> UniqRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub fn new(data: Vec<Range<T>>) -> Self {
        let ranges = Ranges::<T>::new(data);

        UniqRanges {
            ranges
        }
    }

    /// Make the UniqRanges<T> consistent
    /// 
    /// # Info 
    /// 
    /// By construction, the data are sorted so that it is possible (see the new
    /// method definition above) to merge the overlapping ranges.
    pub fn make_consistent(mut self) -> Self {
        self.ranges = self.ranges.make_consistent();
        self
    }

    pub fn to_nested(self) -> NestedRanges<T> {
        let nested_data = UniqToNestedIter::new(self.ranges)
            .collect::<Vec<_>>();
        NestedRanges::<T>::new(nested_data).make_consistent()
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

impl<T> From<Ranges<T>> for UniqRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn from(ranges: Ranges<T>) -> Self {
        UniqRanges::<T> {
            ranges
        }
    }
}

impl<T> From<UniqRanges<T>> for Ranges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn from(uniq_ranges: UniqRanges<T>) -> Self {
        uniq_ranges.ranges
    }
}

use ndarray::Array2;

impl From<UniqRanges<u64>> for Array2<u64> {
    fn from(input: UniqRanges<u64>) -> Self {
        input.ranges.into()
    }
}

use ndarray::Array1;
fn uniq_ranges_to_array1d<T>(input: UniqRanges<T>) -> Array1<T> 
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::iter::Step
    + std::fmt::Debug {
    let ranges = input.ranges;

    let mut result: Vec<T> = Vec::<T>::new();

    // Add the UNIQs contained into the spatial MOC
    for range in ranges.iter() {
        for uniq_range in range.start..range.end {
            result.push(uniq_range);
        }
    }
    // Get an Array1 from the Vec<i64> without copying any data
    Array1::from_vec(result).to_owned()
}
impl From<UniqRanges<u64>> for Array1<u64> {
    fn from(input: UniqRanges<u64>) -> Self {
        uniq_ranges_to_array1d(input)
    }
}
impl From<UniqRanges<i64>> for Array1<i64> {
    fn from(input: UniqRanges<i64>) -> Self {
        uniq_ranges_to_array1d(input)
    }
}
