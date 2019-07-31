#![feature(test)]
#![feature(step_trait)]
extern crate num;
extern crate rand;
extern crate test;
extern crate rayon;
extern crate ndarray;
extern crate healpix;

pub mod bounded;
mod utils;

mod ranges;
pub mod nestedranges;
pub mod uniqranges;

pub mod nestedranges2d;

use ndarray::Array2;
pub struct RangesPy<T>
{
    pub data: Array2<T>,
}

use crate::nestedranges::NestedRanges;
use crate::uniqranges::UniqRanges;
use num::{Integer, PrimInt};
use bounded::Bounded;

impl<T> RangesPy<T>
where T: Integer + PrimInt
    + Bounded<T>
    + Send + Sync
    + std::fmt::Debug {

    /// Convert the python input ranges to a NestedRanges<T> consistent type.
    pub fn to_nested_ranges(self) -> NestedRanges<T> {
        let data = utils::array2_to_vec_ranges(self.data);

        NestedRanges::<T>::new(data)
    }

    /// Convert the python input ranges to a UniqRanges<T> consistent type.
    pub fn to_uniq_ranges(self) -> UniqRanges<T> {
        let data = utils::array2_to_vec_ranges(self.data);

        UniqRanges::<T>::new(data)
    }
}

#[cfg(test)]
mod tests {
    use crate::nestedranges::NestedRanges;
    use crate::uniqranges::UniqRanges;

    use num::PrimInt;

    #[test]
    fn test_uniq_iter() {
        let simple_nested = NestedRanges::<u64>::new(vec![0..1]);
        let complex_nested = NestedRanges::<u64>::new(vec![7..76]);
        let empty_nested = NestedRanges::<u64>::new(vec![]);

        let simple_uniq = UniqRanges::<u64>::new(vec![4*4.pow(29)..(4*4.pow(29) + 1)]);
        let complex_uniq = UniqRanges::<u64>::new(
            vec![
                (1 + 4*4.pow(27))..(4 + 4*4.pow(27)),
                (2 + 4*4.pow(28))..(4 + 4*4.pow(28)),
                (16 + 4*4.pow(28))..(19 + 4*4.pow(28)),
                (7 + 4*4.pow(29))..(8 + 4*4.pow(29))
            ]
        ).make_consistent();
        let empty_uniq = UniqRanges::<u64>::new(vec![]).make_consistent();

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
        
        let ranges = UniqRanges::<u64>::new(input.clone()).make_consistent();
        let expected = UniqRanges::<u64>::new(input).make_consistent();

        assert_eq!(ranges.to_nested().to_uniq(), expected);
    }
}