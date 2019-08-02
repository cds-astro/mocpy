use num::{Integer, PrimInt};

use crate::ranges::{Ranges, NestedToUniqIter, DepthPixIter};
use crate::uniqranges::UniqRanges;
use crate::bounded::Bounded;

use std::ops::Range;
use std::slice::Iter;

#[derive(Debug, Clone)]
pub struct NestedRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    ranges: Ranges<T>
}

impl<T> NestedRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    pub fn new(data: Vec<Range<T>>) -> Self {
        let ranges = Ranges::<T>::new(data);

        NestedRanges {
            ranges
        }
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
    pub fn divide(mut self, min_depth: i8) -> Self {
        self.ranges = self.ranges.divide(min_depth);
        self
    }

    pub fn to_uniq(self) -> UniqRanges<T> {
        let uniq_data = NestedToUniqIter::new(self.ranges)
            .collect::<Vec<_>>();
        UniqRanges::<T>::new(uniq_data).make_consistent()
    }

    pub fn depth(&self) -> i8 {
        self.ranges.depth()
    }

    pub fn degrade(&mut self, depth: i8) {
        self.ranges.degrade(depth);
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

    pub fn iter_depth_pix(self) -> DepthPixIter<T> {
        DepthPixIter::<T>::new(self.ranges)
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

impl<T> From<Ranges<T>> for NestedRanges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn from(ranges: Ranges<T>) -> Self {
        NestedRanges::<T> {
            ranges
        }
    }
}

impl<T> From<NestedRanges<T>> for Ranges<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {
    fn from(nested_ranges: NestedRanges<T>) -> Self {
        nested_ranges.ranges
    }
}

use ndarray::Array2;

impl From<NestedRanges<u64>> for Array2<u64> {
    fn from(input: NestedRanges<u64>) -> Self {
        input.ranges.into()
    }
}
impl From<NestedRanges<i64>> for Array2<i64> {
    fn from(input: NestedRanges<i64>) -> Self {
        input.ranges.into()
    }
}