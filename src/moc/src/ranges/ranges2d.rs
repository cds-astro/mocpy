
use std::cmp;
use std::ops::Range;
use std::collections::HashSet;

use rayon::prelude::*;

use crate::idx::Idx;
use crate::ranges::{Ranges, SNORanges};

/// Generic operations on a set of Sorted and Non-Overlapping ranges.
/// SNO = Sorted Non-Overlapping
pub trait SNORanges2D<'a, T: Idx, S: Idx>: Sized {

    // TODO: remove from here and put directly in a constructor
    fn make_consistent(self) -> Self;

    fn is_empty(&self) -> bool;

    /// Checks whether a (time, position) tuple lies into the `NestedRanges2D<T, S>`
    ///
    /// # Arguments
    ///
    /// * ``time`` - The time at which the coordinate has been observed
    /// * ``range`` - The spatial pixel of the nested range
    fn contains(&self, time: T, range: &Range<S>) -> bool;

    fn union(&self, other: &Self) -> Self;
    fn intersection(&self, other: &Self) -> Self;
    fn difference(&self, other: &Self) -> Self;
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



#[derive(Debug)]
pub struct Ranges2D<T: Idx, S: Idx> {
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


#[derive(Eq, PartialEq, Debug)]
struct BoundRange<T: Idx> {
    x: T,
    y_idx: usize,
    start: bool
}

impl<T: Idx> BoundRange<T> {
    fn new(x: T, y_idx: usize, start: bool) -> BoundRange<T> {
        BoundRange {
            x,
            y_idx,
            start
        }
    }
}

use std::cmp::Ordering;
impl<T: Idx> Ord for BoundRange<T> {
    fn cmp(&self, other: &BoundRange<T>) -> Ordering {
        // Notice that the we flip the ordering on costs.
        // In case of a tie we compare positions - this step is necessary
        // to make implementations of `PartialEq` and `Ord` consistent.
        self.x.cmp(&other.x)
          .then_with(|| self.start.cmp(&other.start))
    }
}

// `PartialOrd` needs to be implemented as well.
impl<T: Idx> PartialOrd for BoundRange<T> {
    fn partial_cmp(&self, other: &BoundRange<T>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: Idx, S: Idx> Ranges2D<T, S> {
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

    /// Private method for merging touching time intervals
    /// that share a common spatial coverages
    fn compress(&mut self, time_ranges: Vec<Range<T>>, spatial_coverages: Vec<Ranges<S>>) {
        // Vector storing the final time ranges
        let mut res_t = Vec::<Range<T>>::with_capacity(time_ranges.len());
        // Vector storing the final spatial coverages
        let mut res_s = Vec::<Ranges<S>>::with_capacity(spatial_coverages.len());

        // The first tuple is kept and we will start iterating
        // from the second tuple.
        let mut prev_s = spatial_coverages.first().unwrap().clone();
        let mut prev_t = time_ranges.first().unwrap().clone();
        // At this point:
        // * There is at least 2 (time, space) tuples.
        // * Time ranges are not overlapping because (see the step 2. for more details)
        //   They can touch to each other.
        //
        // We loop over the time ranges to merge the touching ones which have
        // equal spatial coverages.
        for (cur_t, cur_s) in time_ranges.into_iter().zip(spatial_coverages.into_iter()).skip(1) {
            match cur_t.start.cmp(&prev_t.end) {
                Ordering::Greater => {
                    // a. Time ranges are not overlapping
                    res_t.push(prev_t);
                    res_s.push(prev_s);

                    prev_t = cur_t;
                    prev_s = cur_s;
                },
                Ordering::Equal => {
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
                },
                Ordering::Less => {
                    // c. They cannot overlap each other
                    unreachable!("no overlapping time ranges at this point");
                }
            }
        }
        res_t.push(prev_t);
        res_s.push(prev_s);

        res_s.shrink_to_fit();
        res_t.shrink_to_fit();

        self.x = res_t;
        self.y = res_s;
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

        let mut i = 0_usize;
        let mut j = 0_usize;

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
            } else if let Some(cur_ranges) = s {
                if !cur_ranges.is_empty() {
                    last_t = Some(c);

                    s_ranges.push(cur_ranges);
                    prev_s = s_ranges.last();
                }
            }
        }

        Ranges2D {
            x: t_ranges,
            y: s_ranges,
        }
    }
}

impl<'a, T: Idx, S: Idx> SNORanges2D<'a, T, S> for Ranges2D<T, S> {

    fn make_consistent(mut self) -> Self {
        if !self.is_empty() {
            let mut sorted_time_bound_ranges = self.x.par_iter()
              .enumerate()
              .map(|(idx, t)| {
                  let start_b = BoundRange::new(t.start, idx, true);
                  let end_b = BoundRange::new(t.end, idx, false);

                  vec![start_b, end_b]
              })
              .flatten()
              .collect::<Vec<_>>();

            (&mut sorted_time_bound_ranges).par_sort_unstable_by(|l, r| {
                l.cmp(&r)
            });

            let mut time_ranges = vec![];
            let mut spatial_coverages = vec![];

            let mut ranges_idx = HashSet::new();
            ranges_idx.insert(0);
            let mut prev_time_bound = sorted_time_bound_ranges
              .first()
              .unwrap()
              .x;
            for time_bound in sorted_time_bound_ranges.iter().skip(1) {
                let cur_time_bound = time_bound.x;

                if cur_time_bound > prev_time_bound && !ranges_idx.is_empty() {
                    let spatial_coverage = ranges_idx.par_iter()
                      .map(|coverage_idx| self.y[*coverage_idx].clone())
                      .reduce(
                          Ranges::<S>::default,
                          |s1, s2| {
                              s1.union(&s2)
                          }
                      );

                    time_ranges.push(prev_time_bound..cur_time_bound);
                    spatial_coverages.push(spatial_coverage);
                }

                if time_bound.start {
                    // A new time range begins, we add its S-MOC idx
                    // to the set
                    ranges_idx.insert(time_bound.y_idx);
                } else {
                    ranges_idx.remove(&time_bound.y_idx);
                }

                prev_time_bound = cur_time_bound;
            }

            // Time ranges are not overlapping anymore but can be next
            // to each other.
            // One must check if they are equal and if so, merge them.
            self.compress(time_ranges, spatial_coverages);
        }

        self
    }

    fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    fn contains(&self, time: T, range: &Range<S>) -> bool {
        // Check whether the time lies in the ranges of the `T` dimension
        let in_first_dim = self.x
          .par_iter()
          .enumerate()
          .filter_map(|(idx, r)| {
              let in_time_range = time >= r.start && time <= r.end;
              if in_time_range {
                  Some(idx)
              } else {
                  None
              }
          })
          .collect::<Vec<_>>();

        if in_first_dim.is_empty() {
            // The time is not contained in any ranges so we simply
            // return false here
            false
        } else if in_first_dim.len() == 1 {
            let idx_first_dim = in_first_dim.first().unwrap();
            // Check whether the pixel coordinate lies in the `S` dimension
            // coverage where the time lies.
            let s_coverage = &self.y[*idx_first_dim];
            s_coverage.contains(range)
        } else {
            unreachable!();
        }
    }

    fn union(&self, other: &Self) -> Self {
        self.merge(other, Self::op_union)
    }

    fn intersection(&self, other: &Self) -> Self {
        self.merge(other, Self::op_intersection)
    }

    fn difference(&self, other: &Self) -> Self {
        self.merge(other, Self::op_difference)
    }
}



impl<T: Idx, S: Idx> PartialEq for Ranges2D<T, S> {
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
    use std::ops::Range;
    use crate::ranges::{Ranges, Idx};
    use crate::ranges::ranges2d::{Ranges2D, SNORanges2D};

    fn creating_ranges<T: Idx, S: Idx>(
        ranges_t: Vec<Range<T>>,
        ranges_s: Vec<Vec<Range<S>>>,
    ) -> Ranges2D<T, S> {
        let mut vec_ranges_s = Vec::<Ranges<S>>::with_capacity(ranges_t.len());

        for range_s in ranges_s.into_iter() {
            vec_ranges_s.push(Ranges::<S>::new_from(range_s));
        }

        Ranges2D::new(ranges_t, vec_ranges_s).make_consistent()
    }

    #[test]
    fn merge_overlapping_ranges() {
        let t = vec![0..15, 0..15, 15..30, 30..45, 15..30];
        let s = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
          .make_consistent();

        let t_expect = vec![0..15, 15..45];
        let s_expect = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }

    // Overlapping time ranges configuration:
    // xxxxxxxxxxx
    // xxxx-------
    #[test]
    fn remove_different_length_time_ranges() {
        let t = vec![0..7, 0..30];
        let s = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
          .make_consistent();

        let t_expect = vec![0..7, 7..30];
        let s_expect = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
            Ranges::<u64>::new_unchecked(vec![16..21])
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }

    // Overlapping time ranges configuration:
    // xxxxxxxxxxx
    // ----xxxx---
    #[test]
    fn remove_different_length_time_ranges2() {
        let t = vec![0..30, 2..10];
        let s = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
          .make_consistent();

        let t_expect = vec![0..2, 2..10, 10..30];
        let s_expect = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18])
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }

    // Overlapping time ranges configuration:
    // xxxxxxx----
    // ----xxxxxxx
    #[test]
    fn remove_different_length_time_ranges3() {
        let t = vec![0..5, 2..10];
        let s = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
          .make_consistent();

        let t_expect = vec![0..2, 2..5, 5..10];
        let s_expect = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
            Ranges::<u64>::new_unchecked(vec![16..21])
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }

    // Overlapping time ranges configuration:
    // xxxxxxxxxxx
    // ----xxxxxxx
    #[test]
    fn remove_different_length_time_ranges4() {
        let t = vec![0..30, 10..30];
        let s = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
          .make_consistent();

        let t_expect = vec![0..10, 10..30];
        let s_expect = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![0..4, 5..21]),
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }
    // No overlapping time ranges
    // xxxxxx----
    // ------xxxx
    #[test]
    fn remove_different_length_time_ranges5() {
        let t = vec![0..5, 5..20];
        let s = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
          .make_consistent();

        let t_expect = vec![0..5, 5..20];
        let s_expect = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![16..21])
        ];
        let coverage_expect = Ranges2D::<u64, u64>::new(t_expect, s_expect);
        assert_eq!(coverage, coverage_expect);
    }

    #[test]
    fn merge_overlapping_ranges_2() {
        let t = vec![0..15, 0..15, 15..30, 30..45, 15..30];
        let s = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![0..4]),
            Ranges::<u64>::new_unchecked(vec![16..21]),
            Ranges::<u64>::new_unchecked(vec![16..21, 25..26]),
        ];
        let coverage = Ranges2D::<u64, u64>::new(t, s)
          .make_consistent();

        let t_expect = vec![0..15, 15..30, 30..45];
        let s_expect = vec![
            Ranges::<u64>::new_unchecked(vec![0..4, 5..16, 17..18]),
            Ranges::<u64>::new_unchecked(vec![0..4, 16..21, 25..26]),
            Ranges::<u64>::new_unchecked(vec![16..21])
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
