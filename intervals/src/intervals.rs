use crate::ranges::Ranges;
use num::{Integer, PrimInt, CheckedAdd, One, Zero};
use crate::bounded::Bounded;

use std::ops::{Range, BitOr};

/// Intervals definition.
/// 
/// Can be either given in the nested or uniq form.
#[derive(Debug)]
pub enum Intervals<T>
where T: Integer + PrimInt + Bounded<T> + std::fmt::Debug  {
    /// Ranges given in nested form
    Nested(Ranges<T>),
    /// Ranges given in uniq form
    Uniq(Ranges<T>),
}

impl<T> From<Vec<Range<T>>> for Intervals<T>
where T: Integer + PrimInt + Bounded<T> + Send + std::fmt::Debug {
    fn from(ranges: Vec<Range<T>>) -> Self {
        Intervals::Nested(Ranges::<T>::new(ranges, None))
    }
}

impl<T> Intervals<T>
where T: Integer + PrimInt + Bounded<T> + BitOr<T> + Send + std::fmt::Debug + 'static {
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
    ///		None
    ///     )
    /// );
    /// // Call `uniq_into_iter` on it to get an iterator
    /// // that will return ranges under the uniq format
    /// // This operation takes the ownership of `intervals_nest`
    /// let uniq_ranges: Vec<_> = intervals_nest
    ///     .uniq_into_iter()
    ///     .collect();
    /// ```
    pub fn uniq_into_iter(self) -> Box<Iterator<Item = Range<T>>> {
        match self {
            Intervals::Nested(ranges) => {
                Box::new(UniqIntervalsIter::new(ranges))
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
    ///     )
    /// );
    /// // Call `nested_into_iter` on it to get an iterator
    /// // that will return ranges under the nested format
    /// // This operation takes the ownership of `intervals_uniq`
    /// let nested_ranges: Vec<_> = intervals_uniq
    ///     .nested_into_iter()
    ///     .collect();
    /// ```
    pub fn nested_into_iter(self) -> Box<Iterator<Item = Range<T>>> {
        match self {
            Intervals::Uniq(ranges) => {
                Box::new(NestedIntervalsIter::new(ranges))
            },
            Intervals::Nested(ranges) => {
                Box::new(ranges.0.into_iter())
            },
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
    ///     )
    /// );
    /// // Call `to_nested` on it to get a new Intervals
    /// // that contains the ranges expressed in the nested format
    /// let intervals_nested = intervals_uniq.to_nested();
    /// ```
    pub fn to_nested(self) -> Intervals<T> {
        match self {
            Intervals::Uniq(ranges) => {
                let data_nested: Vec<_> = NestedIntervalsIter::new(ranges).collect();
                Intervals::Nested(Ranges::<T>::new(data_nested, None))
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
    ///		None
    ///     )
    /// );
    /// // Call `to_uniq` on it to get a new Intervals object
    /// // that contains the ranges expressed in the uniq format
    /// let intervals_uniq = intervals_nested.to_uniq();
    /// ```
    pub fn to_uniq(self) -> Intervals<T> {
        match self {
            Intervals::Nested(ranges) => {
                let uniq_data: Vec<_> = UniqIntervalsIter::new(ranges).collect();
                Intervals::Uniq(Ranges::<T>::new(uniq_data, None))
            },
            Intervals::Uniq(_) => {
                self
            }
        }
    }

    /// Get an iterator over the range references
    pub fn iter<'a>(&'a self) -> Box<Iterator<Item = &Range<T>> + 'a> {
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

    fn merge(&mut self, other: Self, op: &Fn(bool, bool) -> bool) {
        match self {
            Intervals::Nested(ref mut ranges) => {
                match other {
                    Intervals::Nested(other_ranges) => {
                        ranges.merge(other_ranges, op);
                    },
                    Intervals::Uniq(other_uniq_ranges) => {
                        // Convert uniq intervals to nested ones
                        let other_nested_ranges: Vec<_> = NestedIntervalsIter::new(other_uniq_ranges).collect();
                        ranges.merge(Ranges::<T>::new(other_nested_ranges, None), op);
                    }
                }
            },
            Intervals::Uniq(ref mut ranges) => {
                match other {
                    Intervals::Nested(other_nested_ranges) => {
                        // Convert nested intervals to uniq ones
                        let other_uniq_ranges: Vec<_> = UniqIntervalsIter::new(other_nested_ranges).collect();
                        ranges.merge(Ranges::<T>::new(other_uniq_ranges, None), op);
                    },
                    Intervals::Uniq(other_ranges) => {
                        ranges.merge(other_ranges, op);
                    }
                }
            }
        }
    }

    /// Perform the union between self and another Intervals object
    /// 
    /// # Arguments
    /// 
    /// * `other` - Take the ownership of the intervals
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// 
    /// let mut a = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![12..17, 3..5, 5..7, 6..8],
    ///		None
    ///     )
    /// );
    /// let b = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![0..1, 2..5],
    ///		None
    ///     )
    /// );
    /// a.union(b);
    /// ```
    pub fn union(&mut self, other: Self) {
        self.merge(other, &|a, b| a || b);
    }
    
    /// Perform the difference between self and another Intervals object
    /// 
    /// # Arguments
    /// 
    /// * `other` - Take the ownership of the intervals
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// 
    /// let mut a = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![12..17, 3..5, 5..7, 6..8],
    ///		None
    ///     )
    /// );
    /// let b = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![0..1, 2..5],
    ///		None,
    ///     )
    /// );
    /// a.difference(b);
    /// ```
    pub fn difference(&mut self, other: Self) {
        self.merge(other, &|a, b| a && !b);
    }
    
    /// Perform the intersection between self and another Intervals object
    /// 
    /// # Arguments
    /// 
    /// * `other` - Take the ownership of the intervals
    /// 
    /// Example 
    /// 
    /// ```
    /// use intervals::intervals::Intervals;
    /// use intervals::ranges::Ranges;
    /// 
    /// let mut a = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![12..17, 3..5, 5..7, 6..8],
    ///		None,
    ///     )
    /// );
    /// let b = Intervals::Nested(
    ///     Ranges::<u64>::new(
    ///         vec![0..1, 2..5],
    ///		None
    ///     )
    /// );
    /// a.intersection(b);
    /// ```
    pub fn intersection(&mut self, other: Self) {
        self.merge(other, &|a, b| a && b)
    }

    pub fn complement(&mut self) {
        match self {
            Intervals::Nested(ref mut ranges) => {
                ranges.complement()
            },
            Intervals::Uniq(_) => {
                panic!("Cannot be able to complement a intervals expressed as uniq.\n
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
                let shift = ((<T>::MAXDEPTH - depth) << 1) as u32;

                let mut offset: T = One::one();
                offset = offset.unsigned_shl(shift) - One::one();

                let mut mask: T = One::one();
                mask = mask.checked_mul(&!offset).unwrap();

                let adda: T = Zero::zero();
                let mut addb: T = One::one();
                addb = addb.checked_mul(&offset).unwrap();

                let capacity = ranges.0.len();
                let mut result = Vec::<Range<T>>::with_capacity(capacity);

                for range in ranges.0.iter() {
                    let a: T = range.start.checked_add(&adda).unwrap() & mask;
                    let b: T = range.end.checked_add(&addb).unwrap() & mask;

                    if b > a {
                        result.push(a..b);
                    }
                }

                ranges.0 = result;
            },
            Intervals::Uniq(_) => {
                panic!("Can only degrade nested expressed intervals!");
            }
        }
    }
}

pub struct UniqIntervalsIter<T>
where T: Integer + PrimInt + CheckedAdd + Bounded<T> {
    ranges: Ranges<T>,

    depth: i8,
    shift: u32,
    offset: T,
    depth_offset: T,
}

impl<T> UniqIntervalsIter<T>
where  T: Integer + PrimInt + CheckedAdd + Bounded<T> {
    fn new(ranges: Ranges<T>) -> UniqIntervalsIter<T> {
        let depth = 0;
        let shift = ((T::MAXDEPTH - depth) << 1) as u32;

        let mut offset: T = One::one();
        offset = offset.unsigned_shl(shift) - One::one();

        let mut depth_offset: T = One::one();
        depth_offset = depth_offset.unsigned_shl((2 * depth + 2) as u32);

        UniqIntervalsIter {
            ranges,
            
            depth,
            shift,
            offset,
            depth_offset,
        }
    }
}

impl<T> Iterator for UniqIntervalsIter<T>
where T: Integer + PrimInt + CheckedAdd + Bounded<T> + std::fmt::Debug + Send {
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.ranges.is_empty() {
            for range in self.ranges.iter() {
                let t1 = range.start + self.offset;
                let t2 = range.end;
                
                let pix1 = t1.unsigned_shr(self.shift);
                let pix2 = t2.unsigned_shr(self.shift);
                
                let c1 = pix1.unsigned_shl(self.shift);
                let c2 = pix2.unsigned_shl(self.shift);

                if c2 > c1 {
                    self.ranges.difference(
                        Ranges::<T>::new(vec![c1..c2], None)
                    );

                    let e1 = self.depth_offset.checked_add(&pix1).unwrap();
                    let e2 = self.depth_offset.checked_add(&pix2).unwrap();

                    return Some(e1..e2);
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

pub struct NestedIntervalsIter<T>
where T: Integer + PrimInt + CheckedAdd + Bounded<T> {
    ranges: Ranges<T>,
    u: T,
    
    range_index: usize,
}

impl<T> NestedIntervalsIter<T>
where  T: Integer + PrimInt + CheckedAdd + Bounded<T> {
    fn new(ranges: Ranges<T>) -> NestedIntervalsIter<T> {
        let range_index = 0;

        let u = if ranges.0.len() > 0 {
            ranges[range_index].start
        } else {
            Zero::zero()
        };
        NestedIntervalsIter {
            ranges,
            u,
            range_index,
        }
    }
}

impl<T> Iterator for NestedIntervalsIter<T>
where T: Integer + PrimInt + CheckedAdd + Bounded<T> + std::fmt::Debug {
    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.range_index < self.ranges.0.len() {
            let (depth, ipix) = T::pix_depth(self.u);

            let shift = (T::MAXDEPTH as u32 - depth) << 1;

            let one: T = One::one();
            let e1 = ipix
                .unsigned_shl(shift);
            let e2 = ipix
                .checked_add(&one)
                .unwrap()
                .unsigned_shl(shift);

            self.u = self.u
                .checked_add(&one)
                .unwrap();

            let u_end = self.ranges[self.range_index].end;
            if self.u == u_end {
                self.range_index += 1;

                if self.range_index < self.ranges.0.len() {
                    self.u = self.ranges[self.range_index].start;
                }
            }

            return Some(e1..e2)
        }
        None
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
                    let ra_uniq = Ranges::<u64>::new(left.uniq_into_iter().collect(), None);
                    let rb_uniq = Ranges::<u64>::new(right.uniq_into_iter().collect(), None);

                    assert_eq!(ra_uniq.0, rb_uniq.0);
                }
            },
            (Intervals::Uniq(ref ra), Intervals::Uniq(ref rb)) => {
                if !compare_nested {
                    assert_eq!(ra.0, rb.0);
                } else {
                    let ra_nested = Ranges::<u64>::new(left.nested_into_iter().collect(), None);
                    let rb_nested = Ranges::<u64>::new(right.nested_into_iter().collect(), None);

                    assert_eq!(ra_nested.0, rb_nested.0);
                }
            },
            (Intervals::Nested(ref ra), Intervals::Uniq(ref rb)) => {
                if compare_nested {
                    let rb_nested = Ranges::<u64>::new(right.nested_into_iter().collect(), None);
                    assert_eq!(ra.0, rb_nested.0);
                } else {
                    let ra_uniq = Ranges::<u64>::new(left.uniq_into_iter().collect(), None);
                    assert_eq!(ra_uniq.0, rb.0);
                }
            },
            (Intervals::Uniq(ref ra), Intervals::Nested(ref rb)) => {
                if compare_nested {
                    let ra_nested = Ranges::<u64>::new(left.nested_into_iter().collect(), None);
                    assert_eq!(ra_nested.0, rb.0);
                } else {
                    let rb_uniq = Ranges::<u64>::new(right.uniq_into_iter().collect(), None);
                    assert_eq!(ra.0, rb_uniq.0);
                }
            }
        }
    }

    #[test]
    fn merge_range() {
        fn assert_merge(a: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let intervals = Intervals::Nested(Ranges::<u64>::new(a, None));
            let expected_intervals = Intervals::Nested(Ranges::<u64>::new(expected, None));

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
            Ranges::<u64>::new(vec![0..(1<<58)], Some(1))
        );
	let expected_ranges = vec![0..(1<<56), (1<<56)..(1<<57), (1<<57)..3*(1<<56), 3*(1<<56)..(1<<58)];

	if let Intervals::Nested(ref ranges) = &intervals_d0 {
	    assert_eq!(ranges.0, expected_ranges);
	}
    }

    #[test]
    fn test_uniq_iter() {
        let intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![0..1], None)
        );
        let complex_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![7..76], None)
        );
        let empty_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![], None)
        );

        let expected_intervals = Intervals::Uniq(Ranges::<u64>::new(vec![4*4.pow(29)..(4*4.pow(29) + 1)], None));
        let expected_complex_intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                vec![
                    (1 + 4*4.pow(27))..(4 + 4*4.pow(27)),
                    (2 + 4*4.pow(28))..(4 + 4*4.pow(28)),
                    (16 + 4*4.pow(28))..(19 + 4*4.pow(28)),
                    (7 + 4*4.pow(29))..(8 + 4*4.pow(29))
                ],
		None
            )
        );
        let expected_empty_intervals = Intervals::Uniq(Ranges::<u64>::new(vec![], None));

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
		None
            )
        );
        let expected_intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                input_cloned,
		None
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
		None
            )
        );
        let complex_intervals = Intervals::Uniq(
            Ranges::<u64>::new(
                vec![(1 + 4*4.pow(27))..(4 + 4*4.pow(27)),
                    (2 + 4*4.pow(28))..(4 + 4*4.pow(28)),
                    (16 + 4*4.pow(28))..(19 + 4*4.pow(28)),
                    (7 + 4*4.pow(29))..(8 + 4*4.pow(29))
                ],
		None
            )
        );
        let empty_intervals = Intervals::Uniq(
            Ranges::<u64>::new(vec![], None)
        );

        let expected_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![0..1], None)
        );
        let expected_complex_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![7..76], None)
        );
        let expected_empty_intervals = Intervals::Nested(
            Ranges::<u64>::new(vec![], None)
        );

        intervals_eq(intervals, expected_intervals, true);
        intervals_eq(complex_intervals, expected_complex_intervals, true);
        intervals_eq(empty_intervals, expected_empty_intervals, true);
    }

    #[test]
    fn test_union() {
        fn assert_union(a: Vec<Range<u64>>, b: Vec<Range<u64>>, expected: Vec<Range<u64>>) {
            let mut intervals_a = Intervals::Nested(Ranges::<u64>::new(a, None));
            let intervals_b = Intervals::Nested(Ranges::<u64>::new(b, None));
            
            let expected_union_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None)
            );

            intervals_a.union(intervals_b);
            intervals_eq(intervals_a, expected_union_intervals, true);
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
            let mut intervals_a = Intervals::Nested(Ranges::<u64>::new(a, None));
            let intervals_b = Intervals::Nested(Ranges::<u64>::new(b, None));
            
            let expected_inter_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None)
            );

            intervals_a.intersection(intervals_b);
            intervals_eq(intervals_a, expected_inter_intervals, true);
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
            let mut intervals_a = Intervals::Nested(Ranges::<u64>::new(a, None));
            let intervals_b = Intervals::Nested(Ranges::<u64>::new(b, None));
            
            let expected_diff_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None)
            );

            intervals_a.difference(intervals_b);
            intervals_eq(intervals_a, expected_diff_intervals, true);
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
            let mut intervals = Intervals::Nested(
                Ranges::<u64>::new(input, None)
            );

            let expected_intervals = Intervals::Nested(
                Ranges::<u64>::new(expected, None)
            );

            intervals.complement();
            intervals_eq(intervals, expected_intervals, true);
        }

        fn assert_complement_pow_2(input: Vec<Range<u64>>) {
            let input_cloned = input.clone();
            
            let mut intervals = Intervals::Nested(
                Ranges::<u64>::new(input, None)
            );

            let start_intervals = Intervals::Nested(
                Ranges::<u64>::new(input_cloned, None)
            );

            intervals.complement();
            intervals.complement();

            intervals_eq(intervals, start_intervals, true);
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
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 1)], None)
        );
        assert_eq!(i1.depth(), 0);

        let i2 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 3)], None)
        );
        assert_eq!(i2.depth(), 2);

        let i3 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..3*4.pow(29 - 3)], None)
        );
        assert_eq!(i3.depth(), 3);

        let i4 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..12*4.pow(29)], None)
        );
        assert_eq!(i4.depth(), 0);
    }

    #[test]
    fn test_degrade() {
        let mut i1 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 1)], None)
        );
        i1.degrade(0);
        assert_eq!(i1.depth(), 0);

        let mut i2 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 3)], None)
        );
        i2.degrade(1);
        assert_eq!(i2.depth(), 1);

        let mut i3 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..3*4.pow(29 - 3)], None)
        );
        i3.degrade(1);
        assert_eq!(i3.depth(), 1);

        let mut i4 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..12*4.pow(29)], None)
        );
        i4.degrade(0);
        assert_eq!(i4.depth(), 0);

        let mut i5 = Intervals::Nested(
            Ranges::<u64>::new(vec![0..4*4.pow(29 - 3)], None)
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
