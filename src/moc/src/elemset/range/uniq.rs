use std::ops::Range;
use std::slice::Iter;

use crate::idx::Idx;

use super::HpxRanges;
use super::hpx::UniqToHpxIter;

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

    pub fn into_hpx(self) -> HpxRanges<T> {
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

#[cfg(test)]
mod tests {
    use crate::elemset::range::uniq::HpxUniqRanges;

    #[test]
    fn test_uniq_to_hpx() {
        let uniq = HpxUniqRanges::<u64>::new_unchecked(vec![
            72057594037927937..72057594037927938,
            72057594037927938..72057594037927939,
            72057594037927939..72057594037927940,
            288230376151711746..288230376151711747,
            288230376151711747..288230376151711748,
            288230376151711760..288230376151711761,
            288230376151711761..288230376151711762,
            288230376151711762..288230376151711763,
            1152921504606846983..1152921504606846984
        ]);
        let hpx = uniq.into_hpx();
        println!("hpx: {:?}", hpx);
    }
}
