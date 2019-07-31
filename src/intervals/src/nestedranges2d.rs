
use num::{Integer, PrimInt, One};
use crate::bounded::Bounded;
use crate::ranges::ranges2d::Ranges2D;
use crate::ranges::Ranges;
use crate::nestedranges::NestedRanges;

use rayon::prelude::*;

use std::ops::Range;

#[derive(Debug)]
pub struct NestedRanges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
    ranges: Ranges2D<T, S>,
}

impl<T, S> NestedRanges2D<T, S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {

    /// Create a Quantity/Space 2D coverage:
    ///
    /// - A set of quantity points stored in ``x`` expressed at the maximum depth.
    ///   This quantity axe can be a time one or a redshift one etc...
    ///   This will define the first dimension of the object.
    /// - A set of spatial HEALPix cells stored in ``y`` and given at the depth ``d2``
    ///   This will define the second dimension of the object.
    ///
    /// Quantity could refer to a time but also a proper motion, a redshift, etc...
    /// The resulted 2D coverage will be of depth (``d1``, ``d2``)
    /// 
    /// # Precondition
    /// 
    /// - `d1` and `d2` depths must be valid (within [0, <T>::MAXDEPTH])
    /// - `x` and `y` must have the same size.
    pub fn create_quantity_space_coverage(x: Vec<T>, y: Vec<S>, d1: i8, d2: i8) -> NestedRanges2D<T, S> {
        let s1 = ((<T>::MAXDEPTH - d1) << 1) as u32;
        let mut off1: T = One::one();
        off1 = off1.unsigned_shl(s1) - One::one();

        let mut m1: T = One::one();
        m1 = m1.checked_mul(&!off1).unwrap();

        let x = x.into_par_iter()
                .map(|r| {
                    let a: T = r & m1;
                    let b: T = r.checked_add(&One::one())
                                .unwrap()
                                .checked_add(&off1)
                                .unwrap() & m1;
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
                    // We do not want a min_depth along the 2nd dimension
                    // making sure that the created Ranges<S> is valid.
                    Ranges::<S>::new(vec![a..b])
                })
                .collect::<Vec<_>>();

        let ranges = Ranges2D::<T, S>::new(x, y).make_consistent();

        NestedRanges2D {
            ranges
        }
    }

    /// Create a Space/Quantity 2D coverage:
    ///
    /// - A set of spatial HEALPix cells stored in ``x`` and given at the depth ``d1``
    ///   This will define the first dimension of the object.
    /// - A set of quantity points stored in ``y`` expressed at the maximum depth.
    ///   This quantity axe can be a time one or a redshift one etc...
    ///   This will define the second dimension of the object.
    ///
    /// Quantity could refer to a time but also a proper motion, a redshift, etc...
    /// The resulted 2 dim ranges will be of depth (``d1``, ``d2``)
    /// 
    /// # Precondition
    /// 
    /// - `d1` and `d2` depths must be valid (within [0, <T>::MAXDEPTH])
    /// - `x` and `y` must have the same size.
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
                    // We do not want a min_depth along the 2nd dimension
                    // making sure that the created Ranges<S> is valid.
                    Ranges::<S>::new(vec![a..b])
                })
                .collect::<Vec<_>>();

        let ranges = Ranges2D::<T, S>::new(x, y).make_consistent();

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
                // We do not want a min_depth along the 2nd dimension
                // making sure that the created Ranges<S> is valid.
                || Ranges::<S>::new(vec![]),
                |s1, s2| {
                    s1.union(&s2)
                }
            );

        ranges.into()
    }

    // Returns all the 1st dim ranges for which their 2nd dim set of ranges
    // are contained into ``y``
    //
    // Works in the same manner of ``project_on_first_dim``
    pub fn project_on_first_dim(y: &NestedRanges<S>, coverage: &NestedRanges2D<T, S>) -> NestedRanges<T> {
        let coverage = &coverage.ranges;
        let t_ranges = coverage.x
            .par_iter()
            .zip_eq(coverage.y.par_iter())
            // Filter the time ranges to keep only those
            // that lie into ``x``
            .filter_map(|(t, s)| {
                for r in s.iter() {
                    if !y.contains(r) {
                        return None;
                    }
                }
                // The matching 1st dim ranges matching
                // are cloned. We do not want
                // to consume the Range2D
                Some(t.clone())
            })
            .collect::<Vec<_>>();

        NestedRanges::<T>::new(t_ranges).make_consistent()
    }

    /// Get the maximum depth of the 2D ranges along its first and
    /// second dimensions.
    /// 
    /// If the NestedRanges2D is empty the depth returned is set
    /// to (0, 0)
    pub fn depth(&self) -> (i8, i8) {
        let coverage = &self.ranges;
        let y = coverage.y
            .par_iter()
            // Compute the depths of the Ranges<T>
            .map(|ranges| ranges.depth())
            // Get the max of these depths
            .max()
            // If there are no ranges, the max depth
            // along the second dimension is set to 0
            .unwrap_or_else(|| {
                0
            });

        let x = coverage.x
            .par_iter()
            // Compute de depths of the first dimensional ranges
            .map(|range| {
                let res = range.start | range.end;
                let mut depth: i8 = <T>::MAXDEPTH - (res.trailing_zeros() >> 1) as i8;

                if depth < 0 {
                    depth = 0;
                }
                depth
            })
            // Get the max of these depths
            .max()
            // If there are no ranges, the max depth
            // along the first dimension is set to 0
            .unwrap_or_else(|| {
                0
            });

        (x, y)
    }

    /// Returns the minimum value along the `T` dimension
    pub fn t_min(&self) -> Result<T, &'static str> {
        if self.ranges.is_empty() {
            Err("The coverage is empty")
        } else {
            Ok(self.ranges.x[0].start)
        }
    }

    /// Returns the maximum value along the `T` dimension
    pub fn t_max(&self) -> Result<T, &'static str> {
        if self.ranges.is_empty() {
            Err("The coverage is empty")
        } else {
            Ok(self.ranges.x.last().unwrap().end)
        }
    }
}

use ndarray::Array1;
impl From<&NestedRanges2D<u64, u64>> for Array1<i64> {
    /// Create a Array1<i64> from a NestedRanges2D<u64, u64>
    /// 
    /// This is used when storing a STMOC into a FITS file
    /// 
    /// # Info
    /// 
    /// The output Array1 stores the STMOC under the nested format.
    /// Its memory layout contains each time range followed by the
    /// list of space ranges referred to that time range.
    /// Time ranges are negatives so that one can distinguish them
    /// from space ranges.
    /// 
    /// Content example of an Array1 coming from a FITS file:
    /// int64[] = {-1, -3, 3, 5, 10, 12, 13, 18, -5, -6, 0, 1}
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
            for second_dim_range in s.iter() {
                result.push(second_dim_range.start as i64);
                result.push(second_dim_range.end as i64);
            }
        }

        // Get an Array1 from the Vec<i64> without copying any data
        Array1::from_vec(result).to_owned()
    }
}

use crate::utils;
impl From<Array1<i64>> for NestedRanges2D<u64, u64> {
    /// Create a NestedRanges2D<u64, u64> from a Array1<i64>
    /// 
    /// This is used when loading a STMOC from a FITS file
    /// opened with astropy
    /// 
    /// # Precondition
    /// 
    /// The input Array1 stores the STMOC under the nested format.
    /// Its memory layout contains each time range followed by the
    /// list of space ranges referred to that time range.
    /// Time ranges are negatives so that one can distinguish them
    /// from space ranges.
    /// 
    /// Content example of an Array1 coming from a FITS file:
    /// int64[] = {-1, -3, 3, 5, 10, 12, 13, 18, -5, -6, 0, 1}
    /// 
    /// Coverages coming from FITS file are consistent because they
    /// are stored this way.
    fn from(input: Array1<i64>) -> Self {
        let ranges = if input.is_empty() {
            // If the input array is empty
            // then we return an empty coverage
            Ranges2D::<u64, u64>::new(vec![], vec![])
        } else {
            let mut input = input.into_raw_vec();
            let input = utils::unflatten(&mut input);

            let mut t = Vec::<Range<u64>>::new();
            
            let mut cur_s = Vec::<Range<u64>>::new();
            let mut s = Vec::<Ranges<u64>>::new();
            for r in input.into_iter() {
                if r.start < 0 {
                    // First dim range
                    let t_start = (-r.start) as u64;
                    let t_end = (-r.end) as u64;
                    t.push(t_start..t_end);

                    // Push the second dim MOC if there is ranges in it
                    if !cur_s.is_empty() {
                        // Warning: We suppose all the STMOCS read from FITS files
                        // are already consistent when they have been saved!
                        // That is why we do not check the consistency of the MOCs here!
                        s.push(Ranges::<u64>::new(cur_s.clone()));
                        cur_s.clear();
                    }
                } else {
                    // Second dim range
                    let s_start = r.start as u64;
                    let s_end = r.end as u64;
                    cur_s.push(s_start..s_end);
                }
            }

            // Push the last second dim coverage
            s.push(Ranges::<u64>::new(cur_s));

            assert_eq!(t.len(), s.len());
            // No need to make it consistent because it comes
            // from python
            Ranges2D::<u64, u64>::new(t, s)
        };

        NestedRanges2D {
            ranges
        }
    }
}