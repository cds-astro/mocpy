
use rayon::ThreadPoolBuildError;

pub mod idx;
pub mod qty;
pub mod elem;
pub mod elemset;
pub mod ranges;
pub mod moc;

pub mod mocranges2d;
pub mod hpxranges2d;
pub mod moc2d;

pub mod deser;

pub mod utils;

// from_fits
// from_cone
// from_...

/// Init the number of threads for parallel tasks.
/// Must be called only once!
/// If not called, the default number of threads is the nmber of physical core.
/// See [rayon doc](https://docs.rs/rayon/1.5.1/rayon/struct.ThreadPoolBuilder.html)
pub fn init_par(num_threads: usize) -> Result<(), ThreadPoolBuildError> {
  rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global()
}

#[cfg(test)]
mod tests {
    use num::PrimInt;

    use crate::elemset::range::{
        HpxRanges,
        uniq::HpxUniqRanges
    };

    #[test]
    fn test_uniq_iter() {
        let simple_nested = HpxRanges::<u64>::new_unchecked(vec![0..1]);
        let complex_nested = HpxRanges::<u64>::new_unchecked(vec![7..76]);
        let empty_nested = HpxRanges::<u64>::default();

        let simple_uniq = HpxUniqRanges::<u64>::new_unchecked(vec![4 * 4.pow(29)..(4 * 4.pow(29) + 1)]);
        let complex_uniq = HpxUniqRanges::<u64>::new_from_sorted(vec![
            (1 + 4 * 4.pow(27))..(4 + 4 * 4.pow(27)),
            (2 + 4 * 4.pow(28))..(4 + 4 * 4.pow(28)),
            (16 + 4 * 4.pow(28))..(19 + 4 * 4.pow(28)),
            (7 + 4 * 4.pow(29))..(8 + 4 * 4.pow(29)),
        ]);
        let empty_uniq = HpxUniqRanges::<u64>::new_unchecked(vec![]);

        assert_eq!(simple_nested.clone().into_hpx_uniq(), simple_uniq);
        assert_eq!(complex_nested.clone().into_hpx_uniq(), complex_uniq);
        assert_eq!(empty_nested.clone().into_hpx_uniq(), empty_uniq);

        assert_eq!(simple_uniq.into_hpx(), simple_nested);
        assert_eq!(complex_uniq.into_hpx(), complex_nested);
        assert_eq!(empty_uniq.into_hpx(), empty_nested);
    }

    #[test]
    fn test_uniq_nested_conversion() {
        let input = vec![
            1056..1057,
            1057..1058,
            1083..1084,
            1048539..1048540,
            1048574..1048575,
            1048575..1048576,
        ];

        let ranges = HpxUniqRanges::<u64>::new_from_sorted(input.clone());
        let expected = HpxUniqRanges::<u64>::new_from_sorted(input);

        assert_eq!(ranges.into_hpx().into_hpx_uniq(), expected);
    }
}
