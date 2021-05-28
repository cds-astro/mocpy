use std::ops::Range;
use std::slice::Iter;

use num::One;
use ndarray::Array1;

use crate::ranges::Idx;
use crate::mocranges::HpxRanges;
use crate::hpxranges::UniqToHpxIter;


#[derive(Debug)]
pub struct HpxUniqRanges<T: Idx> {
    ranges: HpxRanges<T>,
}

impl<T: Idx> HpxUniqRanges<T> {
    pub fn new_unchecked(data: Vec<Range<T>>) -> Self {
        let ranges = HpxRanges::<T>::new_unchecked(data);
        HpxUniqRanges { ranges }
    }

    pub fn new_from_sorted(data: Vec<Range<T>>) -> Self {
        let ranges = HpxRanges::<T>::new_from_sorted(data);
        HpxUniqRanges { ranges }
    }

    pub fn to_hpx(self) -> HpxRanges<T> {
        let nested_data = UniqToHpxIter::new(self.ranges).collect::<Vec<_>>();
        HpxRanges::<T>::new_from(nested_data)
    }

    pub fn iter(&self) -> Iter<Range<T>> {
        self.ranges.iter()
    }
}

impl<T: Idx> PartialEq for HpxUniqRanges<T> {
    fn eq(&self, other: &Self) -> bool {
        self.ranges == other.ranges
    }
}

impl<T> Eq for HpxUniqRanges<T> where T: Idx {}

impl<T: Idx> From<HpxRanges<T>> for HpxUniqRanges<T> {
    fn from(ranges: HpxRanges<T>) -> Self {
        HpxUniqRanges::<T> { ranges }
    }
}

impl<T: Idx> From<HpxUniqRanges<T>> for HpxRanges<T> {
    fn from(uniq_ranges: HpxUniqRanges<T>) -> Self {
        uniq_ranges.ranges
    }
}

use ndarray::Array2;

impl From<HpxUniqRanges<u64>> for Array2<u64> {
    fn from(input: HpxUniqRanges<u64>) -> Self {
        input.ranges.into()
    }
}

fn uniq_ranges_to_array1d<T: Idx>(input: HpxUniqRanges<T>) -> Array1<T> {
    let ranges = input.ranges;

    let mut result: Vec<T> = Vec::<T>::new();

    // Add the UNIQs contained into the spatial MOC
    for range in ranges.iter() {
        for uniq_range in num::range_step(range.start, range.end, One::one()) {
            result.push(uniq_range);
        }
    }
    // Get an Array1 from the Vec<i64> without copying any data
    let result: Array1<T> = result.into();
    result.to_owned()
}
impl From<HpxUniqRanges<u64>> for Array1<u64> {
    fn from(input: HpxUniqRanges<u64>) -> Self {
        uniq_ranges_to_array1d(input)
    }
}
impl From<HpxUniqRanges<i64>> for Array1<i64> {
    fn from(input: HpxUniqRanges<i64>) -> Self {
        uniq_ranges_to_array1d(input)
    }
}
