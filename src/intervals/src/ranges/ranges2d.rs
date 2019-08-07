use crate::bounded::Bounded;
use num::{Integer, PrimInt};

use crate::ranges::Ranges;
use std::cmp;
use std::ops::Range;

use rayon::prelude::*;

#[derive(Debug)]
pub struct Ranges2D<T, S>
where
    T: Integer + Sync + std::fmt::Debug,
    S: Integer + PrimInt + Bounded<S> + Sync + Send + std::fmt::Debug,
{
    // First dimension
    pub x: Vec<Range<T>>,
    // Second dimension (usually the spatial one)
    // Consecutive S values do not always refer to
    // neighbours spatial cells. Therefore there is a few chance
    // that time coverage will be the same for consecutive
    // spatial cells. So the spatial cells will not be merged
    // a lot.
    pub y: Vec<Ranges<S>>,
}

type Operation<T, S> = fn(
    // Left 2D ranges operand
    &Ranges2D<T, S>,
    // Right 2D ranges operand
    &Ranges2D<T, S>,
    // On rising edge or in one range
    // of the left operand
    bool,
    // On rising edge or in one range
    // of the right operand
    bool,
    // The index of the current range bound
    // of the left operand
    usize,
    // The index of the current range bound
    // of the right operand
    usize,
) -> Option<Ranges<S>>;

impl<T, S> Ranges2D<T, S>
where
    T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
    S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug,
{
    /// Creates a 2D coverage
    ///
    /// # Arguments
    ///
    /// * `t` - A set of ranges constituing the first dimension. This stores
    ///   usually quantities such as times, redshifts or proper motions.
    /// * `s` - A set of 1D coverage constituing the second dimension. This stores
    ///   usually space informations such as HEALPix cell indices under the nested format.
    ///
    /// # Precondition
    ///
    /// ``t`` and ``s`` must have the same length.
    pub fn new(t: Vec<Range<T>>, s: Vec<Ranges<S>>) -> Ranges2D<T, S> {
        Ranges2D { x: t, y: s }
    }

    pub fn make_consistent(mut self) -> Self {
        fn intersect_range<T: Integer + Clone + Copy>(r1: &Range<T>, r2: &Range<T>) -> Option<Range<T>> {
            if r1.start >= r2.end || r2.start >= r1.end {
                // No intersection
                None
            } else {
                let r3_start = std::cmp::max(r1.start, r2.start);
                let r3_end = std::cmp::min(r1.end, r2.end);

                Some(r3_start..r3_end)
            }
        }

        let size = self.x.len();
        // Consistency has sense only from 2 tuples
        if size > 1 {
            // 1. Tuples are joined and sorted along their first dimension.
            let mut t_s = self.x
                .into_par_iter()
                .zip_eq(self.y.into_par_iter())
                .collect::<Vec<_>>();

            (&mut t_s).par_sort_unstable_by(|l, r| {
                l.0.start.cmp(&r.0.start)
            });

            // 2. Times ranges created by a single time at a specific depth
            //    are all of the same length.
            //    Moreover there cannot be time ranges with different starting time
            //    that overlap. In other words, overlapping time ranges must be doublons.
            //    
            //    This step removes the time ranges doublons by doing the union of their
            //    spatial coverages.
            let mut cur_t = &t_s.first().unwrap().0;
            let mut t_idx_to_merge = Vec::<Vec<usize>>::new();
            t_idx_to_merge.push(vec![0]);

            for (i, (t, _)) in t_s.iter().enumerate().skip(1) {
                if t != cur_t {
                    cur_t = t;
                    t_idx_to_merge.push(vec![i]);
                } else {
                    t_idx_to_merge.last_mut()
                        .unwrap()
                        .push(i);
                }
            }

            let (t, s): (Vec<_>, Vec<_>) = t_idx_to_merge
                .into_iter()
                .map(|t_idx| {
                    let first_idx = t_idx.first().unwrap().clone();
                    let t = t_s[first_idx].0.clone();

                    let s = t_idx.into_par_iter()
                        .map(|i| {
                            let r = t_s[i].1.clone();
                            r
                        })
                        .reduce(
                            || Ranges::<S>::new(vec![]),
                            |s1, s2| {
                                s1.union(&s2)
                            }
                        );

                    (t, s)
                })
                .unzip();

            // Vector storing the final time ranges
            let mut res_t = Vec::<Range<T>>::with_capacity(t.len());
            // Vector storing the final spatial coverages
            let mut res_s = Vec::<Ranges<S>>::with_capacity(s.len());

            // The first tuple is kept and we will start iterating
            // from the second tuple.
            let mut prev_s = s[0].clone();
            let mut prev_t = t[0].clone();
            // At this point: 
            // * There is at least 2 (time, space) tuples.
            // * Time ranges are not overlapping because (see the step 2. for more details)
            //   They can touch to each other.
            //
            // We loop over the time ranges to merge the touching ones which have
            // equal spatial coverages.
            for (cur_t, cur_s) in t.into_iter().zip(s.into_iter()).skip(1) {
                if cur_t.start > prev_t.end {
                    // a. Time ranges are not overlapping
                    res_t.push(prev_t);
                    res_s.push(prev_s);

                    prev_t = cur_t;
                    prev_s = cur_s;
                } else if cur_t.start == prev_t.end {
                    // b. Time ranges are touching to each other
                    if cur_s == prev_s {
                        // We merge time intervals if their
                        // spatial coverages are equal
                        prev_t.end = cur_t.end;
                    } else {
                        // If not we push the previous range
                        res_t.push(prev_t);
                        res_s.push(prev_s);

                        prev_t = cur_t;
                        prev_s = cur_s;
                    }
                } else {
                    // c. They cannot overlap each other
                    unreachable!("no overlapping time ranges at this point");
                }
            }
            res_t.push(prev_t);
            res_s.push(prev_s);

            res_s.shrink_to_fit();
            res_t.shrink_to_fit();

            self.x = res_t;
            self.y = res_s;
        }

        self
    }

    fn merge(&self, other: &Self, op: Operation<T, S>) -> Ranges2D<T, S> {
        // Get the ranges from the first dimensions of self and
        // cast them to flat vectors
        let t1 = &self.x;
        let t1_l = t1.len() << 1;

        // Get the first dimension ranges from other
        let t2 = &other.x;
        let t2_l = t2.len() << 1;

        let mut t_ranges = Vec::with_capacity(3 * cmp::max(t1_l, t2_l));
        let mut s_ranges = Vec::with_capacity(3 * cmp::max(t1_l, t2_l));

        let mut i = 0 as usize;
        let mut j = 0 as usize;

        // We will just need a reference to the previous
        // S Ranges because we will only compare it
        // to the current S Ranges.
        // If it is equal, then we do not need to change
        // anything. If not, we have to push the previous
        // S Ranges to the resulting S Ranges stack and set
        // its value to the current S Ranges.
        let mut prev_s: Option<&Ranges<S>> = None;

        let mut last_t = None;

        while i < t1_l || j < t2_l {
            let (c, s) = if i == t1_l {
                let v2 = if j & 0x1 != 0 {
                    t2[j >> 1].end
                } else {
                    t2[j >> 1].start
                };

                let c = v2;

                let in_t1 = false;
                let on_rising_edge_t2 = (j & 0x1) == 0;
                let in_t2 = on_rising_edge_t2;

                let s = op(self, other, in_t1, in_t2, i, j);

                j += 1;

                (c, s)
            } else if j == t2_l {
                let v1 = if i & 0x1 != 0 {
                    t1[i >> 1].end
                } else {
                    t1[i >> 1].start
                };

                let c = v1;

                let on_rising_edge_t1 = (i & 0x1) == 0;
                let in_t1 = on_rising_edge_t1;
                let in_t2 = false;

                let s = op(self, other, in_t1, in_t2, i, j);

                i += 1;

                (c, s)
            } else {
                let v1 = if i & 0x1 != 0 {
                    t1[i >> 1].end
                } else {
                    t1[i >> 1].start
                };
                let v2 = if j & 0x1 != 0 {
                    t2[j >> 1].end
                } else {
                    t2[j >> 1].start
                };

                let c = cmp::min(v1, v2);

                let on_rising_edge_t1 = (i & 0x1) == 0;
                let on_rising_edge_t2 = (j & 0x1) == 0;
                let in_t1 = (on_rising_edge_t1 && c == v1) | (!on_rising_edge_t1 && c < v1);
                let in_t2 = (on_rising_edge_t2 && c == v2) | (!on_rising_edge_t2 && c < v2);

                let s = op(self, other, in_t1, in_t2, i, j);

                if c == v1 {
                    i += 1;
                }
                if c == v2 {
                    j += 1;
                }

                (c, s)
            };

            if let Some(prev_ranges) = prev_s {
                if let Some(cur_ranges) = s {
                    if cur_ranges.is_empty() {
                        // Case of the difference
                        // The difference of two 1D coverage
                        // can be NULL therefore the time range must not
                        // be pushed
                        t_ranges.push(last_t.unwrap()..c);
                        last_t = None;
                        prev_s = None;
                    } else if !prev_ranges.eq(&cur_ranges) {
                        t_ranges.push(last_t.unwrap()..c);
                        last_t = Some(c);

                        s_ranges.push(cur_ranges);
                        prev_s = s_ranges.last();
                    }
                } else {
                    t_ranges.push(last_t.unwrap()..c);
                    last_t = None;
                    prev_s = None;
                }
            } else {
                if let Some(cur_ranges) = s {
                    if !cur_ranges.is_empty() {
                        last_t = Some(c);

                        s_ranges.push(cur_ranges);
                        prev_s = s_ranges.last();
                    }
                }
            }
        }

        Ranges2D {
            x: t_ranges,
            y: s_ranges,
        }
    }

    fn op_union(
        &self,
        other: &Self,
        in_t1: bool,
        in_t2: bool,
        i: usize,
        j: usize,
    ) -> Option<Ranges<S>> {
        if in_t1 && in_t2 {
            let s1 = &self.y[i >> 1];
            let s2 = &other.y[j >> 1];
            Some(s1.union(s2))
        } else if !in_t1 && in_t2 {
            let s2 = &other.y[j >> 1];
            Some(s2.clone())
        } else if in_t1 && !in_t2 {
            let s1 = &self.y[i >> 1];
            Some(s1.clone())
        } else {
            None
        }
    }

    pub fn union(&self, other: &Self) -> Self {
        self.merge(other, Self::op_union)
    }

    fn op_intersection(
        &self,
        other: &Self,
        in_t1: bool,
        in_t2: bool,
        i: usize,
        j: usize,
    ) -> Option<Ranges<S>> {
        if in_t1 && in_t2 {
            let s1 = &self.y[i >> 1];
            let s2 = &other.y[j >> 1];
            Some(s1.intersection(s2))
        } else {
            None
        }
    }

    pub fn intersection(&self, other: &Self) -> Self {
        self.merge(other, Self::op_intersection)
    }

    fn op_difference(
        &self,
        other: &Self,
        in_t1: bool,
        in_t2: bool,
        i: usize,
        j: usize,
    ) -> Option<Ranges<S>> {
        if in_t1 && in_t2 {
            let s1 = &self.y[i >> 1];
            let s2 = &other.y[j >> 1];
            Some(s1.difference(s2))
        } else if in_t1 && !in_t2 {
            let s1 = &self.y[i >> 1];
            Some(s1.clone())
        } else {
            None
        }
    }

    pub fn difference(&self, other: &Self) -> Self {
        self.merge(other, Self::op_difference)
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

impl<T, S> PartialEq for Ranges2D<T, S>
where
    T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
    S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug,
{
    fn eq(&self, other: &Self) -> bool {
        // It is fast to check if the two ranges
        // have the same number of ranges towards their
        // first dimension.
        if self.x.len() != other.x.len() {
            false
        } else {
            // If they have we can check if their 1st dim ranges
            // are equal
            if self.x != other.x {
                false
            } else {
                // In the last step we must verify that
                // each 2nd dim ranges are equal
                for (s1, s2) in self.y.iter().zip(other.y.iter()) {
                    if s1 != s2 {
                        return false;
                    }
                }
                true
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bounded::Bounded;
    use crate::ranges::ranges2d::Ranges2D;
    use crate::ranges::Ranges;

    use num::{Integer, PrimInt};
    use std::ops::Range;

    fn creating_ranges<T, S>(
        ranges_t: Vec<Range<T>>,
        ranges_s: Vec<Vec<Range<S>>>,
    ) -> Ranges2D<T, S>
    where
        T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
        S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug,
    {
        let mut vec_ranges_s = Vec::<Ranges<S>>::with_capacity(ranges_t.len());

        for range_s in ranges_s.into_iter() {
            vec_ranges_s.push(Ranges::<S>::new(range_s).make_consistent());
        }

        Ranges2D::new(ranges_t, vec_ranges_s).make_consistent()
    }

    #[test]
    fn merge_overlapping_ranges() {
        let t = vec![0..15, 0..15, 15..30, 30..45, 15..30];
        let s = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new(vec![16..21]),
            Ranges::<u64>::new(vec![16..21]),
            Ranges::<u64>::new(vec![16..21]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
            .make_consistent();

        let t_expect = vec![0..15, 15..45];
        let s_expect = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new(vec![16..21]),
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }

    #[test]
    fn merge_overlapping_ranges_2() {
        let t = vec![0..15, 0..15, 15..30, 30..45, 15..30];
        let s = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new(vec![0..4]),
            Ranges::<u64>::new(vec![16..21]),
            Ranges::<u64>::new(vec![16..21, 25..26]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
            .make_consistent();

        let t_expect = vec![0..15, 15..30, 30..45];
        let s_expect = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new(vec![0..4, 16..21, 25..26]),
            Ranges::<u64>::new(vec![16..21])
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }

    #[test]
    fn union_ranges_1_3() {
        let a = creating_ranges::<u64, u64>(vec![0..10], vec![vec![16..21]]);
        let b = creating_ranges::<u64, u64>(vec![10..20], vec![vec![16..21]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(vec![0..20], vec![vec![16..21]]);
        assert_eq!(res, c);
    }
    #[test]
    fn union_ranges_1_3_bis() {
        let a = creating_ranges::<u64, u64>(vec![0..10], vec![vec![16..21]]);
        let b = creating_ranges::<u64, u64>(vec![10..20], vec![vec![16..22]]);

        let c = a.union(&b);

        let res =
            creating_ranges::<u64, u64>(vec![0..10, 10..20], vec![vec![16..21], vec![16..22]]);
        assert_eq!(res, c);
    }
    #[test]
    fn union_ranges_covering() {
        let a = creating_ranges::<u64, u64>(vec![0..10], vec![vec![16..21]]);
        let b = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(
            vec![0..9, 9..10, 10..20],
            vec![vec![16..21], vec![0..21], vec![0..17]],
        );
        assert_eq!(res, c);
    }

    #[test]
    fn empty_range_union() {
        let a = creating_ranges::<u64, u64>(vec![0..1], vec![vec![42..43]]);
        let b = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(vec![0..1, 9..20], vec![vec![42..43], vec![0..17]]);
        assert_eq!(res, c);
    }

    #[test]
    fn empty_range_union_bis() {
        let b = creating_ranges::<u64, u64>(vec![0..9], vec![vec![0..20]]);
        let a = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(vec![0..9, 9..20], vec![vec![0..20], vec![0..17]]);
        assert_eq!(res, c);
    }

    #[test]
    fn complex_union() {
        let a = creating_ranges::<u64, u64>(
            vec![0..2, 3..5, 8..9, 13..14],
            vec![vec![2..3], vec![2..3], vec![5..6], vec![7..8]],
        );
        let b = creating_ranges::<u64, u64>(
            vec![1..4, 6..7, 9..10, 11..12],
            vec![vec![0..3], vec![5..6], vec![5..6], vec![10..13]],
        );

        let result = a.union(&b);
        let expected = creating_ranges::<u64, u64>(
            vec![0..1, 1..4, 4..5, 6..7, 8..10, 11..12, 13..14],
            vec![
                vec![2..3],
                vec![0..3],
                vec![2..3],
                vec![5..6],
                vec![5..6],
                vec![10..13],
                vec![7..8],
            ],
        );

        assert_eq!(expected, result);
    }

    #[test]
    fn test_intersection() {
        let empty = creating_ranges::<u64, u64>(vec![], vec![]);

        let a = creating_ranges::<u64, u64>(vec![1..4, 6..7], vec![vec![0..3], vec![5..10]]);

        let b = creating_ranges::<u64, u64>(vec![2..3, 6..7], vec![vec![0..5], vec![7..11]]);

        let a_inter_b = a.intersection(&b);
        let expect_a_inter_b =
            creating_ranges::<u64, u64>(vec![2..3, 6..7], vec![vec![0..3], vec![7..10]]);

        assert_eq!(expect_a_inter_b, a_inter_b);

        let b_inter_empty = b.intersection(&empty);
        let expect_b_inter_empty = creating_ranges::<u64, u64>(vec![], vec![]);

        assert_eq!(b_inter_empty, expect_b_inter_empty);

        let empty_inter_a = empty.intersection(&a);
        let expect_empty_inter_a = creating_ranges::<u64, u64>(vec![], vec![]);

        assert_eq!(empty_inter_a, expect_empty_inter_a);

        let empty_inter_empty = empty.intersection(&empty);
        let expect_empty_inter_empty = creating_ranges::<u64, u64>(vec![], vec![]);

        assert_eq!(empty_inter_empty, expect_empty_inter_empty);
    }

    #[test]
    fn test_difference() {
        let empty = creating_ranges::<u64, u64>(vec![], vec![]);

        let a = creating_ranges::<u64, u64>(vec![1..4, 6..7], vec![vec![0..3], vec![5..10]]);

        let b = creating_ranges::<u64, u64>(vec![2..3, 6..7], vec![vec![0..5], vec![7..11]]);

        let c = creating_ranges::<u64, u64>(vec![0..3, 3..7], vec![vec![0..7], vec![0..12]]);

        let a_diff_b = a.difference(&b);
        let expect_a_diff_b = creating_ranges::<u64, u64>(
            vec![1..2, 3..4, 6..7],
            vec![vec![0..3], vec![0..3], vec![5..7]],
        );

        assert_eq!(expect_a_diff_b, a_diff_b);

        let b_diff_empty = b.difference(&empty);

        assert_eq!(b_diff_empty, b);

        let empty_diff_a = empty.difference(&a);
        let expect_empty_diff_a = creating_ranges::<u64, u64>(vec![], vec![]);

        assert_eq!(empty_diff_a, expect_empty_diff_a);

        let empty_diff_empty = empty.difference(&empty);
        let expect_empty_diff_empty = creating_ranges::<u64, u64>(vec![], vec![]);

        assert_eq!(empty_diff_empty, expect_empty_diff_empty);

        let b_diff_c = b.difference(&c);
        let expect_b_diff_c = creating_ranges::<u64, u64>(vec![], vec![]);

        assert_eq!(b_diff_c, expect_b_diff_c);
    }
}
