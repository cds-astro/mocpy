use crate::ranges::Ranges;
use crate::ranges::{UniqToNestedIter,
                    NestedToUniqIter,
                    DepthPixIter};
use num::{Integer, PrimInt, Zero};
use crate::bounded::Bounded;

use std::ops::{Range, BitOr};

/// Intervals definition.
/// 
/// Can be either given in the nested or uniq form.
#[derive(Debug)]
pub enum Intervals<T>
where T: Integer + PrimInt + Bounded<T> + Sync + std::fmt::Debug + Send {
    /// Ranges given in nested form
    Nested(Ranges<T>),
    /// Ranges given in uniq form
    Uniq(Ranges<T>),
}

impl<T> Intervals<T>
where T: Integer + PrimInt + Bounded<T> + BitOr<T> + Send + Sync + std::fmt::Debug + 'static {
    /// Get an iterator returning uniq ranges
    /// 
    /// # Arguments
    /// 
    /// * `self` - Take the ownership of the intervals
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// // Create a set of ranges under the nested format
    /// let intervals_nest = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![0..2],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// // Call `uniq_into_iter` on it to get an iterator
    /// // that will return ranges under the uniq format
    /// // This operation takes the ownership of `intervals_nest`
    /// let uniq_ranges: Vec<_> = intervals_nest
    ///     .uniq_into_iter()
    ///     .collect();
    /// ```
    pub fn uniq_into_iter(self) -> Box<dyn Iterator<Item = Range<T>>> {
        match self {
            Intervals::Nested(ranges) => {
                Box::new(NestedToUniqIter::new(ranges))
            },
            Intervals::Uniq(ranges) => {
                Box::new(ranges.0.into_iter())
            },
        }
    }

    /// Get an iterator returning nested ranges
    /// 
    /// # Arguments
    /// 
    /// * `self` - Take the ownership of the intervals
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// // Create a set of ranges under the uniq format
    /// let intervals_uniq = Intervals::Uniq(
    ///     Ranges::<u64>::new(
    ///         vec![(4*4_u64.pow(12) + 100)..(4*4_u64.pow(12) + 150)],
    ///         None,
    ///         true
    ///     )
    /// );
    /// // Call `nested_into_iter` on it to get an iterator
    /// // that will return ranges under the nested format
    /// // This operation takes the ownership of `intervals_uniq`
    /// let nested_ranges: Vec<_> = intervals_uniq
    ///     .nested_into_iter()
    ///     .collect();
    /// ```
    pub fn nested_into_iter(self) -> Box<dyn Iterator<Item = Range<T>>> {
        match self {
            Intervals::Uniq(ranges) => {
                Box::new(UniqToNestedIter::new(ranges))
            },
            Intervals::Nested(ranges) => {
                Box::new(ranges.0.into_iter())
            },
        }
    }

    pub fn depthpix_into_iter(self) -> Box<dyn Iterator<Item = (i8, T)>> {
        if let Intervals::Nested(ranges) = self {
            Box::new(DepthPixIter::new(ranges))
        } else {
            unreachable!()
        }
    }

    /// Convert the intervals to nested intervals
    /// 
    /// # Arguments
    /// 
    /// * `self` - Take the ownership of the intervals
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// // Create a set of ranges under the uniq format
    /// let intervals_uniq = Intervals::Uniq(
    ///     Ranges::<u64>::new(
    ///         vec![(4*4_u64.pow(12) + 100)..(4*4_u64.pow(12) + 150)],
    ///         None,
    ///         true,
    ///     )
    /// );
    /// // Call `to_nested` on it to get a new Intervals
    /// // that contains the ranges expressed in the nested format
    /// let intervals_nested = intervals_uniq.to_nested();
    /// ```
    pub fn to_nested(self) -> Intervals<T> {
        match self {
            Intervals::Uniq(ranges) => {
                let data_nested: Vec<_> = UniqToNestedIter::new(ranges).collect();
                Intervals::Nested(Ranges::<T>::new(data_nested, None, true))
            },
            Intervals::Nested(_) => {
                self
            }
        }
    }

    /// Convert the intervals to uniq intervals
    /// 
    /// # Arguments
    /// 
    /// * `self` - Take the ownership of the intervals
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// // Create a set of ranges under the nested format
    /// let intervals_nested = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![80..1520, 2..3],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// // Call `to_uniq` on it to get a new Intervals object
    /// // that contains the ranges expressed in the uniq format
    /// let intervals_uniq = intervals_nested.to_uniq();
    /// ```
    pub fn to_uniq(self) -> Intervals<T> {
        match self {
            Intervals::Nested(ranges) => {
                let uniq_data: Vec<_> = NestedToUniqIter::new(ranges).collect();
                Intervals::Uniq(Ranges::<T>::new(uniq_data, None, true))
            },
            Intervals::Uniq(_) => {
                self
            }
        }
    }

    /// Get an iterator over the range references
    pub fn iter<'a>(&'a self) -> Box<dyn Iterator<Item = &Range<T>> + 'a> {
        match self {
            Intervals::Nested(ref ranges) => {
                Box::new(ranges.0.iter())
            },
            Intervals::Uniq(ref ranges) => {
                Box::new(ranges.0.iter())
            },
        }
    }

    /// Check whether the Intervals object is empty (i.e. contains no ranges)
    pub fn is_empty(&self) -> bool {
        match self {
            Intervals::Nested(ref ranges) | Intervals::Uniq(ref ranges) => {
                ranges.is_empty()
            },
        }
    }

    fn merge(&self, other: &Self, op: impl Fn(bool, bool) -> bool) -> Result<Self, &'static str> {
        match self {
            Intervals::Nested(ref ranges) => {
                match other {
                    Intervals::Nested(ref other_ranges) => {
                        let result = ranges.merge(other_ranges, op);
                        Ok(Intervals::<T>::Nested(result))
                    },
                    Intervals::Uniq(_) => {
                        Err("self uses the nested scheme whereas \
                             other uses the uniq scheme")
                    }
                }
            },
            Intervals::Uniq(ref ranges) => {
                match other {
                    Intervals::Nested(_) => {
                        Err("self uses the uniq scheme whereas \
                             other uses the nested scheme")
                    },
                    Intervals::Uniq(other_ranges) => {
                        let result = ranges.merge(other_ranges, op);
                        Ok(Intervals::<T>::Uniq(result))
                    }
                }
            }
        }
    }

    /// Perform the union between self and another Intervals object
    /// 
    /// # Arguments
    /// 
    /// * `other` - Another intervals to perform the union with.
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// 
    /// let a = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![12..17, 3..5, 5..7, 6..8],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// let b = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![0..1, 2..5],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// let result = a.union(&b);
    /// ```
    pub fn union(&self, other: &Self) -> Self {
        self.merge(other, |a, b| a || b).unwrap()
    }
    
    /// Perform the difference between self and another Intervals object
    /// 
    /// # Arguments
    /// 
    /// * `other` - Another intervals to perform the intersection with.
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// 
    /// let a = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![12..17, 3..5, 5..7, 6..8],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// let b = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![0..1, 2..5],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// let result = a.difference(&b);
    /// ```
    pub fn difference(&self, other: &Self) -> Self {
        self.merge(&other, |a, b| a && !b).unwrap()
    }
    
    /// Perform the intersection between self and another Intervals object
    /// 
    /// # Arguments
    /// 
    /// * `other` - Another intervals to perform the difference with.
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// 
    /// let a = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![12..17, 3..5, 5..7, 6..8],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// let b = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![0..1, 2..5],
    ///		    None,
    ///         true
    ///     )
    /// );
    /// let result = a.intersection(&b);
    /// ```
    pub fn intersection(&self, other: &Self) -> Self {
        self.merge(&other, |a, b| a && b).unwrap()
    }

    pub fn complement(&self) -> Result<Self, &'static str> {
        match self {
            Intervals::Nested(ref ranges) => {
                let ranges_stern = ranges.complement();
                let result = Intervals::<T>::Nested(ranges_stern);
                Ok(result)
            },
            Intervals::Uniq(_) => {
                Err("Cannot be able to complement a intervals expressed as uniq. \
                     Please convert it to nested intervals.")
            },
        }
    }

    pub fn depth(&self) -> i8 {
        let total: T = self.iter().fold(Zero::zero(), |res, x| {
            res | x.start | x.end
        });

        let mut depth: i8 = <T>::MAXDEPTH - (total.trailing_zeros() >> 1) as i8;

        if depth < 0 {
            depth = 0;
        }
        depth
    }

    pub fn degrade(&mut self, depth: i8) {
        match self {
            Intervals::Nested(ref mut ranges) => {
                ranges.degrade(depth);
            },
            Intervals::Uniq(_) => {
                panic!("Can only degrade nested expressed intervals!");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bounded::Bounded;
    use crate::intervals::Intervals;
    use crate::ranges::Ranges;

    use num::PrimInt;
    use std::ops::Range;
    use rand::Rng;

    fn intervals_eq(left: Intervals<u64>, right: Intervals<u64>, compare_nested: bool) {
        match (&left, &right) {
            (Intervals::Nested(ref ra), Intervals::Nested(ref rb)) => {
                if compare_nested {
                    assert_eq!(ra.0, rb.0);
                } else {
                    let ra_uniq = Ranges::<u64>::new(left.uniq_into_iter().collect(), None, true);
                    let rb_uniq = Ranges::<u64>::new(right.uniq_into_iter().collect(), None, true);

                    assert_eq!(ra_uniq.0, rb_uniq.0);
                }
            },
            (Intervals::Uniq(ref ra), Intervals::Uniq(ref rb)) => {
                if !compare_nested {
                    assert_eq!(ra.0, rb.0);
                } else {
                    let ra_nested = Ranges::<u64>::new(left.nested_into_iter().collect(), None, true);
                    let rb_nested = Ranges::<u64>::new(right.nested_into_iter().collect(), None, true);

                    assert_eq!(ra_nested.0, rb_nested.0);
                }
            },
            (Intervals::Nested(ref ra), Intervals::Uniq(ref rb)) => {
                if compare_nested {
                    let rb_nested = Ranges::<u64>::new(right.nested_into_iter().collect(), None, true);
                    assert_eq!(ra.0, rb_nested.0);
                } else {
                    let ra_uniq = Ranges::<u64>::new(left.uniq_into_iter().collect(), None, true);
                    assert_eq!(ra_uniq.0, rb.0);
                }
            },
            (Intervals::Uniq(ref ra), Intervals::Nested(ref rb)) => {
                if compare_nested {
                    let ra_nested = Ranges::<u64>::new(left.nested_into_iter().collect(), None, true);
                    assert_eq!(ra_nested.0, rb.0);
                } else {
                    let rb_uniq = Ranges::<u64>::new(right.uniq_into_iter().collect(), None, true);
                    assert_eq!(ra.0, rb_uniq.0);
                }
            }
        }
    }

    #[test]
    fn merge_range() {
        fn assert_merge(a: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let intervals = Intervals::Nested(Ranges::<u64>::new(a, None, true));
            let expected_intervals = Intervals::Nested(Ranges::<u64>::new(expected, None, true));

            intervals_eq(intervals, expected_intervals, true);
        }

        assert_merge(vec![12..17, 3..5, 5..7, 6..8], vec![3..8, 12..17]);
        assert_merge(vec![0..1, 2..5], vec![0..1, 2..5]);
        assert_merge(vec![], vec![]);
        assert_merge(vec![0..6, 7..9, 8..13], vec![0..6, 7..13]);
    }

    #[test]
    fn merge_range_min_depth() {
	let intervals_d0 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..(1<<58)], Some(1), true)
        );
	let expected_ranges = vec![0..(1<<56), (1<<56)..(1<<57), (1<<57)..3*(1<<56), 3*(1<<56)..(1<<58)];

	if let Intervals::Nested(ref ranges) = &intervals_d0 {
	    assert_eq!(ranges.0, expected_ranges);
	}
    }

    #[test]
    fn test_uniq_iter() {
        let intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![0..1], None, true)
        );
        let complex_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![7..76], None, true)
        );
        let empty_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![], None, true)
        );

        let expected_intervals = Intervals::Uniq(Ranges::<u64>::new(vec![4*4.pow(29)..(4*4.pow(29) + 1)], None, true));
        let expected_complex_intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                vec![
                    (1 + 4*4.pow(27))..(4 + 4*4.pow(27)),
                    (2 + 4*4.pow(28))..(4 + 4*4.pow(28)),
                    (16 + 4*4.pow(28))..(19 + 4*4.pow(28)),
                    (7 + 4*4.pow(29))..(8 + 4*4.pow(29))
                ],
		        None,
                true
            )
        );
        let expected_empty_intervals = Intervals::Uniq(Ranges::<u64>::new(vec![], None, true));

        intervals_eq(intervals, expected_intervals, false);
        intervals_eq(complex_intervals, expected_complex_intervals, false);
        intervals_eq(empty_intervals, expected_empty_intervals, false);
    }

    #[test]
    fn test_uniq_nested_conversion() {
        let input = vec![1056..1057, 1057..1058, 1083..1084, 1048539..1048540, 1048574..1048575, 1048575..1048576];
        let input_cloned = input.clone();
        
        let intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                input,
		        None,
                true
            )
        );
        let expected_intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                input_cloned,
		        None,
                true
            )
        );

        let nested_intervals = intervals.to_nested();
        let final_intervals = nested_intervals.to_uniq();

        intervals_eq(expected_intervals, final_intervals, false);
    }

    #[test]
    fn test_nested_iter() {
        let intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                vec![4*4.pow(29)..(4*4.pow(29) + 1)],
		        None,
                true
            )
        );
        let complex_intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                vec![(1 + 4*4.pow(27))..(4 + 4*4.pow(27)),
                    (2 + 4*4.pow(28))..(4 + 4*4.pow(28)),
                    (16 + 4*4.pow(28))..(19 + 4*4.pow(28)),
                    (7 + 4*4.pow(29))..(8 + 4*4.pow(29))
                ],
		        None,
                true
            )
        );
        let empty_intervals = Intervals::Uniq(
            Ranges::<u64>::new(vec![], None, true)
        );

        let expected_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![0..1], None, true)
        );
        let expected_complex_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![7..76], None, true)
        );
        let expected_empty_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![], None, true)
        );

        intervals_eq(intervals, expected_intervals, true);
        intervals_eq(complex_intervals, expected_complex_intervals, true);
        intervals_eq(empty_intervals, expected_empty_intervals, true);
    }

    #[test]
    fn test_union() {
        fn assert_union(a: Vec<Range<u64>>, b: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let intervals_a = Intervals::Nested(Ranges::<u64>::new(a, None, true));
            let intervals_b = Intervals::Nested(Ranges::<u64>::new(b, None, true));
            
            let expected_union_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None, true)
            );

            let result = intervals_a.union(&intervals_b);
            intervals_eq(result, expected_union_intervals, true);
        }

        assert_union(vec![12..17, 3..5, 5..7, 6..8], vec![0..1, 2..5], vec![0..1, 2..8, 12..17]);
        assert_union(vec![12..17, 3..5, 5..7, 6..8], vec![12..17, 3..5, 5..7, 6..8], vec![3..8, 12..17]);
        assert_union(vec![], vec![], vec![]);
        assert_union(vec![12..17], vec![], vec![12..17]);
        assert_union(vec![0..1, 2..3, 4..5], vec![1..22], vec![0..22]);
        assert_union(vec![0..10], vec![15..22], vec![0..10, 15..22]);
    }

    #[test]
    fn test_intersection() {
        fn assert_intersection(a: Vec<Range<u64>>, b: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let intervals_a = Intervals::Nested(Ranges::<u64>::new(a, None, true));
            let intervals_b = Intervals::Nested(Ranges::<u64>::new(b, None, true));
            
            let expected_inter_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None, true)
            );

            let result = intervals_a.intersection(&intervals_b);
            intervals_eq(result, expected_inter_intervals, true);
        }

        assert_intersection(vec![12..17, 3..5, 5..7, 6..8], vec![0..1, 2..5], vec![3..5]);
        assert_intersection(vec![], vec![0..1, 2..5], vec![]);
        assert_intersection(vec![], vec![], vec![]);
        assert_intersection(vec![2..6], vec![0..3, 4..8], vec![2..3, 4..6]);
        assert_intersection(vec![2..6], vec![2..6, 7..8], vec![2..6]);
        assert_intersection(vec![10..11], vec![10..11], vec![10..11]);
    }

    #[test]
    fn test_difference() {
        fn assert_difference(a: Vec<Range<u64>>, b: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let intervals_a = Intervals::Nested(Ranges::<u64>::new(a, None, true));
            let intervals_b = Intervals::Nested(Ranges::<u64>::new(b, None, true));
            
            let expected_diff_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None, true)
            );

            let result = intervals_a.difference(&intervals_b);
            intervals_eq(result, expected_diff_intervals, true);
        }

        assert_difference(vec![0..20], vec![5..7], vec![0..5, 7..20]);
        assert_difference(vec![0..20], vec![0..20], vec![]);
        assert_difference(vec![0..20], vec![], vec![0..20]);
        assert_difference(vec![0..20], vec![19..22], vec![0..19]);
        assert_difference(vec![0..20], vec![25..27], vec![0..20]);
        assert_difference(vec![0..20], vec![1..2, 3..4, 5..6], vec![0..1, 2..3, 4..5, 6..20]);
    }

    #[test]
    fn test_complement() {
        fn assert_complement(input: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let intervals = Intervals::Nested(
                Ranges::<u64>::new(input, None, true)
            );

            let expected_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None, true)
            );

            let result = intervals.complement().unwrap();
            intervals_eq(result, expected_intervals, true);
        }

        fn assert_complement_pow_2(input: Vec<Range<u64>>) {
            let input_cloned = input.clone();
            
            let intervals = Intervals::Nested(
                Ranges::<u64>::new(input, None, true)
            );

            let start_intervals = Intervals::Nested(
                Ranges::<u64>::new(input_cloned, None, true)
            );

            let result = intervals.complement().unwrap();
            let result = result.complement().unwrap();

            intervals_eq(result, start_intervals, true);
        }

        assert_complement(vec![5..10], vec![0..5, 10..u64::MAXPIX]);
        assert_complement_pow_2(vec![5..10]);

        assert_complement(vec![0..10], vec![10..u64::MAXPIX]);
        assert_complement_pow_2(vec![0..10]);

        assert_complement(vec![0..1, 2..3, 4..5, 6..u64::MAXPIX], vec![1..2, 3..4, 5..6]);
        assert_complement_pow_2(vec![0..1, 2..3, 4..5, 6..u64::MAXPIX]);

        assert_complement(vec![], vec![0..u64::MAXPIX]);
        assert_complement_pow_2(vec![]);
    }

    #[test]
    fn test_depth() {
        let i1 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 1)], None, true)
        );
        assert_eq!(i1.depth(), 0);

        let i2 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 3)], None, true)
        );
        assert_eq!(i2.depth(), 2);

        let i3 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..3*4.pow(29 - 3)], None, true)
        );
        assert_eq!(i3.depth(), 3);

        let i4 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..12*4.pow(29)], None, true)
        );
        assert_eq!(i4.depth(), 0);
    }

    #[test]
    fn test_degrade() {
        let mut i1 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 1)], None, true)
        );
        i1.degrade(0);
        assert_eq!(i1.depth(), 0);

        let mut i2 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 3)], None, true)
        );
        i2.degrade(1);
        assert_eq!(i2.depth(), 1);

        let mut i3 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..3*4.pow(29 - 3)], None, true)
        );
        i3.degrade(1);
        assert_eq!(i3.depth(), 1);

        let mut i4 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..12*4.pow(29)], None, true)
        );
        i4.degrade(0);
        assert_eq!(i4.depth(), 0);

        let mut i5 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 3)], None, true)
        );
        i5.degrade(5);
        assert_eq!(i5.depth(), 2);
    }

    #[test]
    fn test_uniq_decompose() {
        macro_rules! uniq_to_pix_depth {
            ($t:ty, $size:expr) => {
                let mut rng = rand::thread_rng();

                (0..$size).for_each(|_| {
                    let depth = rng.gen_range(0, <$t>::MAXDEPTH) as u32;

                    let npix = 12 * 4.pow(depth);
                    let pix = rng.gen_range(0, npix);

                    let uniq = 4*4.pow(depth) + pix;
                    assert_eq!(<$t>::pix_depth(uniq), (depth, pix));
                });
            };
        }

        uniq_to_pix_depth!(u128, 10000);
        uniq_to_pix_depth!(u64, 10000);
        uniq_to_pix_depth!(u32, 10000);
        uniq_to_pix_depth!(u8, 10000);
    }

    use test::Bencher;

    #[bench]
    fn bench_uniq_to_depth_pix(b: &mut Bencher) {
        let mut rng = rand::thread_rng();
        let n = test::black_box(100000);    

        let uniq: Vec<u64> = (0..n).map(|_| {
            let depth = rng.gen_range(0, 30);

            let npix = 12 * 4.pow(depth);
            let pix = rng.gen_range(0, npix);

            let u = 4 * 4.pow(depth) + pix;
            u
        }).collect();

        b.iter(|| {
            uniq.iter().fold(0, |a, b| a + (u64::pix_depth(*b).0 as u64))
        });
    }
}
