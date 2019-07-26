use crate::ranges::Ranges;
use crate::NestedRanges;
use num::{One, Integer, PrimInt};
use crate::bounded::Bounded;

use std::ops::Range;
use std::cmp;
use std::mem;

use rayon::prelude::*;

#[derive(Debug)]
pub struct Ranges2D<T, S>
where T: Integer + Clone + Copy + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Clone + Copy + Sync + Send + std::fmt::Debug {
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
    usize)
-> Option<Ranges<S>>;

impl<T, S> Ranges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {

    fn new(mut t: Vec<Range<T>>,
           mut s: Vec<Ranges<S>>,
           min_depth: Option<i8>,
           make_consistent: bool) -> Result<Ranges2D<T, S>, &'static str> {
        
        if t.len() != s.len() {
            return Err("Each 1st dimension range must refer to an \
                        associated set of 2nd dim ranges. \
                        Therefore `s` and `t` must be of the same size.");
        }

        let (t_ranges, s_ranges) = if make_consistent {
            // This method is responsible for creating a 2D ranges object that have the
            // following behaviour:
            // * A range on the first dimension (containing elements of type T) refers to
            //   one and only one set of ranges of the second dimension (containing elements of type S).
            // * So to speak, ranges of the first dimension that overlap are merged only if they refer
            //   to equal set of ranges of the second dimension. Otherwise, they are kept separately.
            // * A range on the first dimension referring to an empty set of ranges is discarded.
            
            // First step, the 1st dim ranges are sorted
            let mut join_ts = t.into_par_iter().zip_eq(s.into_par_iter())
                               .collect::<Vec<_>>();

            (&mut join_ts).par_sort_unstable_by(|left, right| left.0.start.cmp(&right.0.start));

            let mut t_current = &join_ts[0].0;

            let mut s_to_union = Vec::<Vec<usize>>::new();
            s_to_union.push(vec![0]);

            for (i, (t, _)) in join_ts
                .iter()
                .skip(1)
                .enumerate() {
                    if !t.eq(&t_current) {
                        t_current = t;
                        s_to_union.push(vec![i]);
                    } else {
                        s_to_union.last_mut()
                            .unwrap()
                            .push(i);
                    }
                }

            let (t, s): (Vec<_>, Vec<_>) = s_to_union.into_iter()
                .map(|v| {
                    let first_id: usize = v.first().unwrap().clone();
                    let ref t = join_ts[first_id].0;

                    let s = v.into_par_iter()
                            .map(|i| {
                                let r = join_ts[i].1.clone();
                                r
                            })
                            .reduce(
                                || Ranges::<S>::new(vec![], None, false),
                                |mut s1, mut s2| {
                                    (&mut s1).union_mut(&mut s2);
                                    s1
                                }
                            );

                    (t.clone(), s)
                })
                .unzip();

            // These stacks will contain the final 1st dim set of ranges
            // and for each of them, their corresponding 2nd dim set of ranges.
            let mut res_t = Vec::<Range<T>>::with_capacity(t.len());
            let mut res_s = Vec::<Ranges<S>>::with_capacity(s.len());

            // We begin to loop over the 1st dim ranges to find the first
            // range that refer to a non null 2nd dim set of ranges.
            let mut start_idx = 0 as usize;
            // This range is stored in ``prev_t`` along with its corresponding
            // set of 2nd dim ranges ``prev_s``
            let (prev_t, prev_s): (Option<Range<T>>, Option<Ranges<S>>) = {
                let mut first_t: Option<Range<T>> = None;
                let mut first_s: Option<Ranges<S>> = None;
                for (cur_t, cur_s) in t.iter().zip(s.iter()) {
                    start_idx = start_idx + 1;
                    if !cur_s.is_empty() {
                        first_t = Some(cur_t.clone());
                        first_s = Some(cur_s.clone());
                        break;
                    }
                }
                (first_t, first_s)
            };

            // If there is at least one valid 1st dim range (i.e. one that contains
            // a non null set of 2nd dim ranges).
            if prev_s.is_some() {
                let mut prev_s = prev_s.unwrap();
                let mut prev_t = prev_t.unwrap();
                // We continue looping over the 1st dim ranges from the ``start_idx`` index.
                for (cur_t, mut cur_s) in t.into_iter().skip(start_idx)
                                       .zip(s.into_iter().skip(start_idx)) {
                    // We discard ranges that are not valid
                    // (i.e. those associated with an empty set of 2nd dim ranges).
                    if !cur_s.is_empty() {
                        // If ``prev_t`` and ``cur_t`` are touching to each other
                        // and the 2nd dim ranges are equal
                        if cur_t.start == prev_t.end && cur_s == prev_s {
                            prev_t.end = cur_t.end;
                        // If ``prev_t`` and ``cur_t`` are overlapping
                        } else if cur_t.start < prev_t.end {
                            // We merge the ranges only if their 2nd set of ranges
                            // are matching
                            if cur_s == prev_s {
                                prev_t.end = cmp::max(prev_t.end, cur_t.end);
                            } else if cur_t == prev_t {
                                prev_s.union_mut(&mut cur_s);
                            // If they are not we can be in two cases: 
                            // 1. The second range ``cur_t`` is not contained into ``prev_t``
                            // xxxx|----
                            // --|yyyyyy
                            // Therefore we will get:
                            // xx|z|yyyy
                            // where the chain of ``z`` is the union between the first and the second
                            // set of 2nd dim ranges.
                            //
                            // 2. The second range ``cur_t`` is contained into ``prev_t``
                            // xxxxxx|---
                            // -|yy|-----
                            // In this case we need to get:
                            // x|zz|x|---
                            // where the chain of ``z`` is the union between the first and the second
                            // set of 2nd dim ranges.
                            } else {
                                let t1 = prev_t.start..cur_t.start;

                                let e1 = cmp::min(cur_t.end, prev_t.end);
                                let e2 = cmp::max(cur_t.end, prev_t.end);
                                let t2 = cur_t.start..e1;
                                let t3 = e1..e2;
                                
                                if t1.start != t1.end {
                                    res_t.push(t1);
                                    // We push ``x``
                                    res_s.push(prev_s.clone());
                                }

                                res_t.push(t2.clone());
                                // We push ``z`` aka the union between
                                // ``x`` and ``y``
                                let z = prev_s.union(&cur_s);
                                res_s.push(z);

                                if t3.start != t3.end {
                                    res_t.push(t3.clone());
                                    // Depending on whether we lie on the first
                                    // or the second case, we push either ``x``
                                    // or ``y``
                                    if e2 == prev_t.end {
                                        // 2nd case
                                        res_s.push(prev_s.clone());
                                    } else {
                                        // 1st case
                                        res_s.push(cur_s.clone());
                                        prev_s = cur_s;
                                    }
                                    prev_t = t3;
                                } else {
                                    prev_t = t2;
                                }
                            }
                        // If ``prev_t`` and ``cur_t`` are not overlapping or
                        // if they are touching but the 2nd ranges are not equal
                        } else {
                            res_t.push(prev_t);
                            res_s.push(prev_s);

                            prev_t = cur_t;
                            prev_s = cur_s;
                        }
                    }
                }

                res_t.push(prev_t);
                res_s.push(prev_s);
            }

            res_s.shrink_to_fit();
            res_t.shrink_to_fit();

            (res_t, res_s)
        } else {
            (t, s)
        };

        Ok(Ranges2D {
            x: t_ranges,
            y: s_ranges,
        })
    }
    
    fn merge(&self, other: &Self, op: Operation<T, S>) -> Ranges2D<T, S> {
        // Unflatten a stack containing T-typed elements
        // to a stack of Range<T> types elements without
        // copying the data.
        fn unflatten<T>(input: &mut Vec<T>) -> Vec<Range<T>> {
            let mut owned_input = Vec::<T>::new();
            // We swap the content refered by input with a new
            // allocated vector.
            // This fix the problem when ``input`` is freed by reaching out
            // the end of the caller scope.
            std::mem::swap(&mut owned_input, input);

            let len = owned_input.len() >> 1;
            let cap = owned_input.capacity();
            let ptr = owned_input.as_mut_ptr() as *mut Range<T>;
            
            mem::forget(owned_input);

            let result = unsafe {
                Vec::from_raw_parts(ptr, len, cap)
            };
            
            result
        }

        // Flatten a stack containing Range<T> typed elements to a stack containing
        // the start followed by the end value of the set of ranges (i.e. a Vec<T>).
        // This does a copy of the data. This is necessary because we do not want to
        // modify ``self`` as well as ``other`` and we want to return the result of 
        // the union of the two Ranges2D.
        fn flatten<T>(input: &Vec<Range<T>>) -> Vec<T>
        where T: Integer + Clone + Copy {
            input.clone()
                 .into_iter()
                 // Convert Range<T> to Vec<T> containing
                 // the start and the end values of the range.
                 .map(|r| vec![r.start, r.end])
                 // We can call flatten on a iterator containing other
                 // iterators (or collections in our case).
                 .flatten()
                 // Collect to get back a newly created Vec<T> 
                 .collect()
        }

        let sentinel = <T>::MAXPIX + One::one();
        // Get the ranges from the first dimensions of self and
        // cast them to flat vectors
        let mut t1 = flatten(&self.x);
        // Push the sentinel
        t1.push(sentinel);
        // Get the first dimension ranges from other
        let mut t2 = flatten(&other.x);
        // Push the sentinel
        t2.push(sentinel);

        let mut t_ranges: Vec<T> = Vec::<T>::new();
        let mut s_ranges: Vec<Ranges<S>> = Vec::<Ranges<S>>::new();
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

        while i < t1.len() || j < t2.len() {
            let c = cmp::min(t1[i], t2[j]);
            // If the two Ranges2D have been processed
            // then we break the loop
            if c == sentinel {
                break;
            }

            let on_rising_edge_t1 = (i & 0x1) == 0;
            let on_rising_edge_t2 = (j & 0x1) == 0;
            let in_t1 = (on_rising_edge_t1 && c == t1[i]) | (!on_rising_edge_t1 && c < t1[i]);
            let in_t2 = (on_rising_edge_t2 && c == t2[j]) | (!on_rising_edge_t2 && c < t2[j]);

            let s = op(self, other, in_t1, in_t2, i, j);

            if let Some(prev_ranges) = prev_s {
                if let Some(cur_ranges) = s {
                    if !prev_ranges.eq(&cur_ranges) {
                        t_ranges.push(c);
                        t_ranges.push(c);
                        s_ranges.push(cur_ranges);
                        prev_s = s_ranges.last();
                    }
                } else {
                    t_ranges.push(c);
                    prev_s = None;
                }
            } else {
                if let Some(cur_ranges) = s {
                    t_ranges.push(c);
                    s_ranges.push(cur_ranges);
                    prev_s = s_ranges.last();
                }
            }

            if c == t1[i] {
                i += 1;
            }
            if c == t2[j] {
                j += 1;
            }
        }

        let t_ranges = unflatten(&mut t_ranges);
        Ranges2D {
            x: t_ranges,
            y: s_ranges,
        }
    }

    fn op_union(&self, other: &Self, in_t1: bool, in_t2: bool, i: usize, j: usize) -> Option<Ranges<S>> {
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

    fn union(&self, other: &Self) -> Self {
        self.merge(other, Self::op_union)
    }

    fn op_intersection(&self, other: &Self, in_t1: bool, in_t2: bool, i: usize, j: usize) -> Option<Ranges<S>> {
        if in_t1 && in_t2 {
            let s1 = &self.y[i >> 1];
            let s2 = &other.y[j >> 1];
            Some(s1.intersection(s2))
        } else {
            None
        }
    }

    fn intersection(&self, other: &Self) -> Self {
        self.merge(other, Self::op_intersection)
    }

    fn op_difference(&self, other: &Self, in_t1: bool, in_t2: bool, i: usize, j: usize) -> Option<Ranges<S>> {
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

    fn difference(&self, other: &Self) -> Self {
        self.merge(other, Self::op_difference)
    }

    fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

impl<T, S> PartialEq for Ranges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
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
                for (s1, s2) in self.y.iter()
                           .zip(other.y.iter()) {
                    if s1 != s2 {
                        return false;
                    }
                }
                true
            }
        }
    }
}

use ndarray::Array1;
impl From<&NestedRanges2D<u64, u64>> for Array1<i64> {
    fn from(input: &NestedRanges2D<u64, u64>) -> Self {
        let ranges = &input.ranges;

        let first_dim_ranges = &ranges.x;
        let second_dim_ranges = &ranges.y;

        let mut result: Vec<i64> = Vec::<i64>::new();

        // Iterate over the tuples (time range, spatial moc associated)
        for (t, s) in first_dim_ranges.iter().zip(second_dim_ranges.iter()) {
            // 1. Append the time range. The opposite is taken so that one can
            //    recognize it is a first dimensional range
            result.push(-(t.start as i64));
            result.push(-(t.end as i64));

            // 2. Append the spatial ranges describing the spatial coverage
            //    associated to the above time range.
            for second_dim_range in &s.0 {
                result.push(second_dim_range.start as i64);
                result.push(second_dim_range.end as i64);
            }
        }

        // Get an Array1 from the Vec<i64> without copying any data
        Array1::from_vec(result).to_owned()
    }
}

#[derive(Debug)]
pub struct NestedRanges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
    pub ranges: Ranges2D<T, S>,
}

impl<T, S> NestedRanges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {

    // Create a Time/Redshift/..., Space coverage with:
    //
    // - A set of quantity points stored in ``x`` expressed at the maximum depth.
    //   This quantity axe can be a time one or a redshift one etc...
    //   This will define the first dimension of the object.
    // - A set of spatial HEALPix cells stored in ``y`` and given at the depth ``d2``
    //   This will define the second dimension of the object.
    //
    // The resulted 2 dim ranges will be of depth (``d1``, ``d2``)
    pub fn create_quantity_space_coverage(x: Vec<T>, y: Vec<S>, d1: i8, d2: i8) -> NestedRanges2D<T, S> {
        let s1 = ((<T>::MAXDEPTH - d1) << 1) as u32;
        let mut off1: T = One::one();
        off1 = off1.unsigned_shl(s1) - One::one();

        let mut m1: T = One::one();
        m1 = m1.checked_mul(&!off1).unwrap();

        let x = x.into_par_iter()
                .map(|r| {
                    let a: T = r & m1;
                    let b: T = r.checked_add(&off1).unwrap() & m1;
                    a..b
                })
                .collect::<Vec<_>>();

        let s2 = ((<S>::MAXDEPTH - d2) << 1) as u32;
        let y = y.into_par_iter()
                .map(|r| {
                    let a = r.unsigned_shl(s2);
                    let b = r.checked_add(&One::one())
                             .unwrap()
                             .unsigned_shl(s2);
                    Ranges::<S>::new(vec![a..b], None, false)
                })
                .collect::<Vec<_>>();

        let ranges = Ranges2D::<T, S>::new(x, y, None, true).unwrap();
        NestedRanges2D {
            ranges
        }
    }

    // Create a Space, Time/Redshift/... coverage with:
    //
    // - A set of spatial HEALPix cells stored in ``x`` and given at the depth ``d1``
    //   This will define the first dimension of the object.
    // - A set of quantity points stored in ``y`` expressed at the maximum depth.
    //   This quantity axe can be a time one or a redshift one etc...
    //   This will define the second dimension of the object.
    //
    // The resulted 2 dim ranges will be of depth (``d1``, ``d2``)
    pub fn create_space_quantity_coverage(x: Vec<T>, y: Vec<S>, d1: i8, d2: i8) -> NestedRanges2D<T, S> {
        // The spatial dimension as the first one
        let s1 = ((<S>::MAXDEPTH - d1) << 1) as u32;
        let x = x.into_par_iter()
                .map(|r| {
                    let a = r.unsigned_shl(s1);
                    let b = r.checked_add(&One::one())
                             .unwrap()
                             .unsigned_shl(s1);
                    a..b
                })
                .collect::<Vec<_>>();

        // The quantity dimension as the second one
        let s2 = ((<S>::MAXDEPTH - d2) << 1) as u32;
        let mut off2: S = One::one();
        off2 = off2.unsigned_shl(s2) - One::one();

        let mut m2: S = One::one();
        m2 = m2.checked_mul(&!off2).unwrap();

        let y = y.into_par_iter()
                .map(|r| {
                    let a: S = r & m2;
                    let b: S = r.checked_add(&off2).unwrap() & m2;
                    Ranges::<S>::new(vec![a..b], None, false)
                })
                .collect::<Vec<_>>();

        let ranges = Ranges2D::<T, S>::new(x, y, None, true).unwrap();
        NestedRanges2D {
            ranges
        }
    }

    // Returns all the 2nd dim ranges for which their 1st dim range
    // is contained into ``x``
    //
    // This method checks for all the 1st dim ranges that
    // lie into the range set ``x``.
    // It then returns the union of the 2nd dim set of ranges
    // of the matching ranges (i.e. the ones contained into the
    // range set ``x``).
    pub fn project_on_second_dim(x: &NestedRanges<T>, coverage: &NestedRanges2D<T, S>) -> NestedRanges<S> {
        let coverage = &coverage.ranges;
        let ranges = coverage.x
            .par_iter()
            .zip_eq(coverage.y.par_iter())
            // Filter the time ranges to keep only those
            // that lie into ``x``
            .filter_map(|(t, s)| {
                if x.contains(t) {
                    Some(s.clone())
                } else {
                    None
                }
            })
            // Compute the union of all the 2nd dim ranges
            // that have been kept
            .reduce(
                || Ranges::<S>::new(vec![], None, false),
                |s1, s2| {
                    s1.union(&s2)
                }
            );

        NestedRanges {
            ranges
        }
    }

    // Returns all the 1st dim ranges for which their 2nd dim set of ranges
    // are contained into ``y``
    //
    // Works in the same manner of ``project_on_first_dim``
    pub fn project_on_first_dim(y: &NestedRanges<S>, coverage: NestedRanges2D<T, S>) -> NestedRanges<T> {
        let coverage = coverage.ranges;
        let t_ranges = coverage.x
            .into_par_iter()
            .zip_eq(coverage.y.par_iter())
            // Filter the time ranges to keep only those
            // that lie into ``x``
            .filter_map(|(t, s)| {
                for r in s.0.iter() {
                    if !y.contains(r) {
                        return None;
                    }
                }
                // The matching 1st dim ranges matching
                // are cloned. We do not want
                // to consume the Range2D
                Some(t)
            })
            .collect::<Vec<_>>();

        NestedRanges::<T>::new(t_ranges, None, true)
    }

    /// Get the maximum depth of the 2D ranges along its first and
    /// second dimensions.
    pub fn depth(&self) -> (i8, i8) {
        let coverage = &self.ranges;
        let y = coverage.y
            .par_iter()
            .map(|ranges| ranges.depth())
            .max()
            .unwrap_or_else(|| {
                0
            });

        let x = Ranges::<T>::new(coverage.x.clone(), None, false).depth();

        dbg!(x, y)
    }
}

#[cfg(test)]
mod tests {
    use crate::bounded::Bounded;
    use crate::ranges2d::Ranges2D;
    use crate::ranges::Ranges;

    use num::{PrimInt, Integer};
    use std::ops::Range;

    fn creating_ranges<T, S>(ranges_t: Vec<Range<T>>, ranges_s: Vec<Vec<Range<S>>>) -> Ranges2D<T, S>
    where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
          S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
        let mut vec_ranges_s = Vec::<Ranges<S>>::with_capacity(ranges_t.len());

        for range_s in ranges_s.into_iter() {
            vec_ranges_s.push(Ranges::<S>::new(range_s, None, true));
        }

        Ranges2D::new(ranges_t, vec_ranges_s, None, true).unwrap()
    }

    #[test]
    fn merge_overlapping_ranges() {
        let ranges_t = vec![0..15, 7..14, 16..17, 18..19, 19..25];
        let ranges_s = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![16..21], None, true),
            Ranges::<u64>::new(vec![16..21], None, true),
            Ranges::<u64>::new(vec![16..21, 25..26], None, true),
        ];
        Ranges2D::<u64, u64>::new(ranges_t, ranges_s, None, true).unwrap();
    }

    #[test]
    fn merge_overlapping_ranges_with_empty() {
        let ranges_t = vec![0..15, 7..14, 16..17, 18..19, 19..25];
        let ranges_s = vec![
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![0..4, 5..16, 17..18], None, true),
            Ranges::<u64>::new(vec![], None, true),
            Ranges::<u64>::new(vec![16..21], None, true),
            Ranges::<u64>::new(vec![16..21, 25..26], None, true),
        ];
        Ranges2D::<u64, u64>::new(ranges_t, ranges_s, None, true).unwrap();
    }

    #[test]
    fn creating_empty_ranges() {
        let ranges_t = vec![0..15, 7..14];
        let ranges_s = vec![
            Ranges::<u64>::new(vec![], None, true),
            Ranges::<u64>::new(vec![], None, true),
        ];
        let ranges_2d = Ranges2D::<u64, u64>::new(ranges_t, ranges_s, None, true).unwrap();
        assert!(ranges_2d.is_empty());
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

        let res = creating_ranges::<u64, u64>(vec![0..10, 10..20], vec![vec![16..21], vec![16..22]]);
        assert_eq!(res, c);
    }
    #[test]
    fn union_ranges_covering() {
        let a = creating_ranges::<u64, u64>(vec![0..10], vec![vec![16..21]]);
        let b = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(
            vec![0..9, 9..10, 10..20],
            vec![vec![16..21], vec![0..21], vec![0..17]]
        );
        assert_eq!(res, c);
    }
    
    #[test]
    fn empty_range_union() {
        let a = creating_ranges::<u64, u64>(vec![0..1], vec![vec![]]);
        assert!(a.is_empty());
        let b = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(
            vec![9..20],
            vec![vec![0..17]]
        );
        assert_eq!(res, c);
    }

    #[test]
    fn empty_range_union_bis() {
        let b = creating_ranges::<u64, u64>(vec![0..1], vec![vec![]]);
        assert!(b.is_empty());
        let a = creating_ranges::<u64, u64>(vec![9..20], vec![vec![0..17]]);

        let c = a.union(&b);

        let res = creating_ranges::<u64, u64>(
            vec![9..20],
            vec![vec![0..17]]
        );
        assert_eq!(res, c);
    }

    #[test]
    fn complex_union() {
        let a = creating_ranges::<u64, u64>(
            vec![0..2, 3..5, 8..9, 13..14],
            vec![vec![2..3], vec![2..3], vec![5..6], vec![7..8]]
        );
        let b = creating_ranges::<u64, u64>(
            vec![1..4, 6..7, 9..10, 11..12],
            vec![vec![0..3], vec![5..6], vec![5..6], vec![10..13]]
        );

        let result = a.union(&b);
        let expected = creating_ranges::<u64, u64>(
            vec![0..1, 1..4, 4..5, 6..7, 8..10, 11..12, 13..14],
            vec![vec![2..3], vec![0..3], vec![2..3], vec![5..6], vec![5..6], vec![10..13], vec![7..8]]
        );

        assert_eq!(expected, result);
    }
}