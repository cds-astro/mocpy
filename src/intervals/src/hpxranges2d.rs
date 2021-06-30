
use std::slice;
use std::ops::Range;
use std::convert::{From, TryFrom};

use num::One;
use rayon::prelude::*;

use crate::idx::Idx;
use crate::utils;
use crate::moc::{
    ZSorted, NonOverlapping,
    RangeMOCIntoIterator,
    CellMOCIterator, CellMOCIntoIterator,
    CellOrCellRangeMOCIterator,
    range::{RangeMOC, RangeMocIter}
};
use crate::moc2d::{
    HasTwoMaxDepth, MOC2Properties,
    RangeMOC2Iterator,
    range::RangeMOC2Elem
};
use crate::qty::{MocQty, Hpx, Time};
use crate::ranges::{SNORanges, ranges2d::SNORanges2D, Ranges};
use crate::mocranges::{HpxRanges, MocRanges};
use crate::mocranges2d::Moc2DRanges;

/// Declaration of the ST-MOC type
pub type TimeSpaceMoc<T, S> = HpxRanges2D::<T, Time<T>, S>;

// Just to be able to define specific methods on this struct
#[derive(Debug)]
pub struct HpxRanges2D<TT: Idx, T: MocQty<TT>, ST: Idx>(pub Moc2DRanges<TT, T, ST, Hpx<ST>>);

impl<TT: Idx> HpxRanges2D<TT, Time<TT>, TT> {
    pub fn time_space_iter<'a>(&'a self, depth_max_t: u8, depth_max_s: u8) -> TimeSpaceRangesIter<'a, TT> {
        // let (depth_max_t, depth_max_s) = self.compute_min_depth();
        TimeSpaceRangesIter {
            depth_max_t,
            depth_max_s,
            it_t: self.0.ranges2d.x.iter().peekable(),
            it_s: self.0.ranges2d.y.iter().peekable(),
        }
    }
}

impl<TT, T, ST> HpxRanges2D<TT, T, ST>
where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
{
    /// Create a new empty `NestedRanges2D<T, S>`
    pub fn new() -> HpxRanges2D<TT, T, ST> {
        let ranges = Moc2DRanges::new(vec![], vec![]);
        HpxRanges2D(ranges)
    }
    
    /// Create a Quantity/Space 2D coverage
    ///
    /// # Arguments
    ///
    /// * `x` - A set of values expressed that will be converted to
    ///   ranges and degraded at the depth ``d1``.
    ///   This quantity axe may refer to a time (expressed in µs), a redshift etc...
    ///   This will define the first dimension of the coverage.
    /// * `y` - A set of spatial HEALPix cell indices at the depth ``d2``.
    ///   This will define the second dimension of the coverage.
    /// * `d1` - The depth of the coverage along its 1st dimension.
    /// * `d2` - The depth of the coverage along its 2nd dimension.
    ///
    /// The resulted 2D coverage will be of depth (``d1``, ``d2``)
    ///
    /// # Precondition
    ///
    /// - `d1` must be valid (within `[0, <T>::MAXDEPTH]`)
    /// - `d2` must be valid (within `[0, <S>::MAXDEPTH]`)
    /// - `x` and `y` must have the same size.
    pub fn create_from_times_positions(
        x: Vec<TT>,
        y: Vec<ST>,
        d1: u8,
        d2: u8,
    ) -> HpxRanges2D<TT, T, ST> {
        let s1 = T::shift_from_depth_max(d1); // ((Self::<T>::MAX_DEPTH - d1) << 1) as u32;
        let mut off1: TT = One::one();
        off1 = off1.unsigned_shl(s1 as u32) - One::one();

        let mut m1: TT = One::one();
        m1 = m1.checked_mul(&!off1).unwrap();

        let x = x
            .into_par_iter()
            .map(|r| {
                let a: TT = r & m1;
                let b: TT = r
                    .checked_add(&One::one())
                    .unwrap()
                    .checked_add(&off1)
                    .unwrap()
                    & m1;
                a..b
            })
            .collect::<Vec<_>>();

        // More generic: Hpx::<ST>::shift_from_depth_max(d2)
        let s2 = ((Hpx::<ST>::MAX_DEPTH - d2) << 1) as u32;
        let y = y
            .into_par_iter()
            .map(|r| {
                let a = r.unsigned_shl(s2);
                let b = r.checked_add(&One::one()).unwrap().unsigned_shl(s2);
                // We do not want a min_depth along the 2nd dimension
                // making sure that the created Ranges<ST> is valid.
                Ranges::<ST>::new_unchecked(vec![a..b])
            })
            .collect::<Vec<_>>();

        let ranges = Ranges2D::<TT, ST>::new(x, y).make_consistent();

        HpxRanges2D(ranges.into())
    }

    /// Create a Quantity/Space 2D coverage
    ///
    /// # Arguments
    ///
    /// * `x` - A set of quantity ranges that will be degraded to the depth ``d1``.
    ///   This quantity axe may refer to a time (expressed in µs), a redshift etc...
    ///   This will define the first dimension of the coverage.
    /// * `y` - A set of spatial HEALPix cell indices at the depth ``d2``.
    ///   This will define the second dimension of the coverage.
    /// * `d2` - The depth of the coverage along its 2nd dimension.
    ///
    /// The resulted 2D coverage will be of depth (``d1``, ``d2``)
    ///
    /// # Precondition
    ///
    /// - `d2` must be valid (within `[0, <S>::MAXDEPTH]`)
    /// - `x` and `y` must have the same size.
    /// - `x` must contain `[a..b]` ranges where `b > a`.
    pub fn create_from_time_ranges_positions(
        x: Vec<Range<TT>>,
        y: Vec<ST>,
        d1: u8,
        d2: u8,
    ) -> HpxRanges2D<TT, T, ST> {
        let s1 = T::shift_from_depth_max(d1);
        let mut off1: TT = One::one();
        off1 = off1.unsigned_shl(s1 as u32) - One::one();

        let mut m1: TT = One::one();
        m1 = m1.checked_mul(&!off1).unwrap();

        let x = x
            .into_par_iter()
            .filter_map(|r| {
                let a: TT = r.start & m1;
                let b: TT = r.end
                    .checked_add(&off1)
                    .unwrap()
                    & m1;
                if b > a {
                    Some(a..b)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        // More generic: Hpx::<ST>::shift_from_depth_max(d2)
        let s2 = ((Hpx::<ST>::MAX_DEPTH - d2) << 1) as u32;
        let y = y
            .into_par_iter()
            .map(|r| {
                let a = r.unsigned_shl(s2);
                let b = r.checked_add(&One::one()).unwrap().unsigned_shl(s2);
                // We do not want a min_depth along the 2nd dimension
                // making sure that the created Ranges<S> is valid.
                Ranges::<ST>::new_unchecked(vec![a..b])
            })
            .collect::<Vec<_>>();

        let ranges = Moc2DRanges::<TT, T, ST, Hpx<ST>>::new(x, y)
            .make_consistent();

        HpxRanges2D(ranges)
    }

    /// Create a Quantity/Space 2D coverage
    ///
    /// # Arguments
    ///
    /// * `x` - A set of quantity ranges that will be degraded to the depth ``d1``.
    ///   This quantity axe may refer to a time (expressed in µs), a redshift etc...
    ///   This will define the first dimension of the coverage.
    /// * `y` - A set of spatial HEALPix cell indices at the depth ``d2``.
    ///   This will define the second dimension of the coverage.
    /// * `d2` - The depth of the coverage along its 2nd dimension.
    ///
    /// The resulted 2D coverage will be of depth (``d1``, ``d2``)
    ///
    /// # Precondition
    ///
    /// - `d2` must be valid (within `[0, <S>::MAXDEPTH]`)
    /// - `x` and `y` must have the same size.
    /// - `x` must contain `[a..b]` ranges where `b > a`.
    pub fn create_from_time_ranges_spatial_coverage(
        x: Vec<Range<TT>>,
        y: Vec<HpxRanges<ST>>,
        d1: u8,
    ) -> HpxRanges2D<TT, T, ST> {
        let s1 = T::shift_from_depth_max (d1) as u32;
        let mut off1: TT = One::one();
        off1 = off1.unsigned_shl(s1) - One::one();

        let mut m1: TT = One::one();
        m1 = m1.checked_mul(&!off1).unwrap();

        let x = x
            .into_par_iter()
            .filter_map(|r| {
                let a: TT = r.start & m1;
                let b: TT = r.end
                    .checked_add(&off1)
                    .unwrap()
                    & m1;
                if b > a {
                    Some(a..b)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        let y = y
            .into_par_iter()
            .map(|r| r.0)
            .collect::<Vec<_>>();

        let ranges = Moc2DRanges::<TT, T, ST, Hpx<ST>>::new(x, y)
            .make_consistent();

        HpxRanges2D(ranges)
    }

    /// Returns the union of the ranges along the `S` axis for which their
    /// `T` ranges intersect ``x``
    ///
    /// # Arguments
    ///
    /// * ``x``- The set of ranges along the `T` axis.
    /// * ``coverage`` - The input coverage
    ///
    /// # Algorithm
    ///
    /// This method checks for all the `T` axis ranges of ``coverage`` that
    /// lie into the range set ``x``.
    ///
    /// It then performs the union of the `S` axis ranges corresponding to the
    /// matching ranges along the `T` axis.
    pub fn project_on_second_dim(
        x: &MocRanges<TT, T>,
        coverage: &HpxRanges2D<TT, T, ST>,
    ) -> HpxRanges<ST> {
        let coverage = &coverage.0.ranges2d;
        let ranges = coverage.x
            .par_iter()
            .zip_eq(coverage.y.par_iter())
            // Filter the time ranges to keep only those
            // that intersects with ``x``
            .filter_map(|(t, s)| {
                if x.intersects(t) {
                    Some(s.clone())
                } else {
                    None
                }
            })
            // Compute the union of all the 2nd
            // dim ranges that have been kept
            .reduce(
                || Ranges::<ST>::default(),
                |s1, s2| s1.union(&s2),
            );

        ranges.into()
    }

    /// Returns the union of the ranges along the `T` axis for which their
    /// `S` ranges is contained in ``y``
    ///
    /// # Arguments
    ///
    /// * ``y``- The set of ranges along the `S` axis.
    /// * ``coverage`` - The input coverage.
    ///
    /// # Algorithm
    ///
    /// This method checks for all the `S` axis ranges of ``coverage`` that
    /// lie into the range set ``y``.
    ///
    /// It then performs the union of the `T` axis ranges corresponding to the
    /// matching ranges along the `S` axis.
    pub fn project_on_first_dim(
        y: &HpxRanges<ST>,
        coverage: &HpxRanges2D<TT, T, ST>,
    ) -> MocRanges<TT, T> {
        let coverage = &coverage.0.ranges2d;
        let t_ranges = coverage.x.par_iter()
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
        // TODO: debug_assert: check is sorted!!
        MocRanges::<TT, T>::new_from_sorted(t_ranges)
    }

    /// Compute the depth of the coverage
    ///
    /// # Returns
    ///
    /// A tuple containing two values:
    ///
    /// * The maximum depth along the `T` axis
    /// * The maximum depth along the `S` axis
    ///
    /// # Info
    ///
    /// If the `NestedRanges2D<T, S>` is empty, the depth returned
    /// is set to (0, 0)
    pub fn compute_min_depth(&self) -> (u8, u8) {
        self.0.compute_min_depth()
    }

    /// Returns the minimum value along the `T` dimension
    ///
    /// # Errors
    ///
    /// When the `NestedRanges2D<T, S>` is empty.
    pub fn t_min(&self) -> Result<TT, &'static str> {
        if self.0.ranges2d.is_empty() {
            Err("The coverage is empty")
        } else {
            Ok(self.0.ranges2d.x[0].start)
        }
    }

    /// Returns the maximum value along the `T` dimension
    ///
    /// # Errors
    ///
    /// When the `NestedRanges2D<T, S>` is empty.
    pub fn t_max(&self) -> Result<TT, &'static str> {
        if self.0.is_empty() {
            Err("The coverage is empty")
        } else {
            Ok(self.0.ranges2d.x.last().unwrap().end)
        }
    }

    /// Performs the union between two `NestedRanges2D<T, S>`
    ///
    /// # Arguments
    ///
    /// * ``other`` - The other `NestedRanges2D<T, S>` to
    ///   perform the union with.
    pub fn union(&self, other: &Self) -> Self {
        let ranges = self.0.union(&other.0);
        HpxRanges2D(ranges)
    }

    /// Performs the intersection between two `NestedRanges2D<T, S>`
    ///
    /// # Arguments
    ///
    /// * ``other`` - The other `NestedRanges2D<T, S>` to
    ///   perform the intersection with.
    pub fn intersection(&self, other: &Self) -> Self {
        let ranges = self.0.intersection(&other.0);
        HpxRanges2D(ranges)
    }

    /// Performs the difference between two `NestedRanges2D<T, S>`
    ///
    /// # Arguments
    ///
    /// * ``other`` - The other `NestedRanges2D<T, S>` to
    ///   perform the difference with.
    pub fn difference(&self, other: &Self) -> Self {
        let ranges = self.0.difference(&other.0);
        HpxRanges2D(ranges)
    }

    /// Check whether a `NestedRanges2D<T, S>` has data in
    /// a (time, ra, dec) tuple.
    ///
    /// # Arguments
    ///
    /// * ``time`` - The time of the tuple
    /// * ``range`` - The position that has been converted to a nested range
    pub fn contains(&self, time: TT, range: &Range<ST>) -> bool {
        self.0.contains(time, range)
    }

    /// Check whether a `NestedRanges2D<T, S>` is empty
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

use std::iter::Peekable;
use ndarray::Array1;
use crate::ranges::ranges2d::Ranges2D;
use crate::moc2d::cell::CellMoc2Iter;
use crate::moc2d::cellcellrange::CellOrCellRangeMoc2Iter;

impl<T: MocQty<u64>> From<&HpxRanges2D<u64, T, u64>> for Array1<i64> {
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
    fn from(input: &HpxRanges2D<u64, T, u64>) -> Self {
        let ranges = &input.0.ranges2d;

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
        let result: Array1<i64> = result.into();
        result.to_owned()
    }
}

impl<TT, T, ST> PartialEq for HpxRanges2D<TT, T, ST>
where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
{
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl<TT, T, ST> Eq for HpxRanges2D<TT, T, ST>
where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
{
}

impl<T:MocQty<u64>> TryFrom<Array1<i64>> for HpxRanges2D<u64, T, u64> {
    type Error = &'static str;
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
    /// Coverages coming from FITS file should be consistent because they
    /// are stored this way.
    ///
    /// # Errors
    ///
    /// * If the number of time ranges do not match the number of
    ///   spatial coverages.
    fn try_from(input: Array1<i64>) -> Result<Self, Self::Error> {
        let ranges = if input.is_empty() {
            // If the input array is empty
            // then we return an empty coverage
            Moc2DRanges::<u64, T, u64, Hpx<u64>>::new(vec![], vec![])
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
                        s.push(Ranges::<u64>::new_unchecked(cur_s.clone()));
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
            s.push(Ranges::<u64>::new_unchecked(cur_s));

            // Propagate invalid Coverage FITS errors.
            if t.len() != s.len() {
                return Err("Number of time ranges and
                    spatial coverages do not match.");
            }
            // No need to make it consistent because it comes
            // from python
            Moc2DRanges::<u64, T, u64, Hpx<u64>>::new(t, s)
        };

        Ok(HpxRanges2D(ranges))
    }
}



impl<T: MocQty<u64>> TryFrom<Array1<u64>> for HpxRanges2D<u64, T, u64> {
    type Error = &'static str;
    /// Create a NestedRanges2D<u64, u64> from a Array1<u64>
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
    /// Coverages coming from FITS file should be consistent because they
    /// are stored this way.
    ///
    /// # Errors
    ///
    /// * If the number of time ranges do not match the number of
    ///   spatial coverages.
    fn try_from(input: Array1<u64>) -> Result<Self, Self::Error> {
        let ranges = if input.is_empty() {
            // If the input array is empty
            // then we return an empty coverage
            Moc2DRanges::<u64, T, u64, Hpx<u64>>::new(vec![], vec![])
        } else {
            let mask = u64::MSB_MASK;
            let mut input = input.into_raw_vec();
            let input = utils::unflatten(&mut input);

            let mut cur_t = Vec::<Range<u64>>::new();
            let mut t: Vec::<Range<u64>> = Vec::with_capacity(input.len() >> 2);

            let mut cur_s = Vec::<Range<u64>>::new();
            let mut s: Vec::<Ranges<u64>> = Vec::with_capacity(input.len());
            for r in input.into_iter() {
                if r.start & r.end & mask == mask {
                    if !cur_s.is_empty(){
                        // Push previous (tranges, srange) tuple
                        for rt in cur_t.drain(..) {
                            t.push(rt);
                            s.push(Ranges::<u64>::new_unchecked(cur_s.clone()));
                            cur_s.clear();
                        }
                        assert!(cur_t.is_empty());
                        // cur_t.clear(); Not needed since we drain
                    }
                    // First dim range
                    let t_start = r.start & mask;
                    let t_end = r.end & mask;
                    cur_t.push(t_start..t_end);
                } else {
                    // Second dim range
                    let s_start = r.start as u64;
                    let s_end = r.end as u64;
                    cur_s.push(s_start..s_end);
                }
            }

            // Push the last (tranges, srange) tuple
            for rt in cur_t {
                t.push(rt);
                s.push(Ranges::<u64>::new_unchecked(cur_s.clone()));
                cur_s.clear();
            }

            // Propagate invalid Coverage FITS errors.
            if t.len() != s.len() {
                return Err("Number of time ranges and
                    spatial coverages do not match.");
            }
            // No need to make it consistent because it comes
            // from python
            Moc2DRanges::<u64, T, u64, Hpx<u64>>::new(t, s)
        };

        Ok(HpxRanges2D(ranges))
    }
}

// The 3 following From contains code redundancy. We probably should do something!

impl<I: RangeMOC2Iterator<
    u64, Time::<u64>, RangeMocIter<u64, Time::<u64>>,
    u64, Hpx::<u64>, RangeMocIter<u64, Hpx::<u64>>,
    RangeMOC2Elem<u64, Time::<u64>, u64, Hpx::<u64>>>
> From<I> for HpxRanges2D<u64, Time<u64>, u64> {

    fn from(it: I) -> Self {
        let mut t = Vec::<Range<u64>>::new();
        let mut s = Vec::<Ranges<u64>>::new();
        for elem in it {
            let (moc_t, moc_s) = elem.mocs();
            /* Simpler but we want to avoid the copy of the s_moc for the last t_range
            for range_t in moc_t.into_range_moc_iter() {
                t.push(range_t);
                s.push(moc_s.moc_ranges().ranges().clone())
            }*/
            let mut it = moc_t.into_range_moc_iter().peekable();
            while it.peek().is_some() {
                let range_t = it.next().unwrap();
                t.push(range_t);
                s.push(moc_s.moc_ranges().ranges().clone())
            }
            if let Some(range_t) = it.next() {
                t.push(range_t);
                s.push(moc_s.into_moc_ranges().into_ranges())
            }
        }
        HpxRanges2D(Moc2DRanges::<u64, Time<u64>, u64, Hpx<u64>>::new(t, s))
    }
}

impl From<CellOrCellRangeMoc2Iter<u64, Time<u64>, u64, Hpx::<u64>>> for HpxRanges2D<u64, Time<u64>, u64> {

    fn from(it: CellOrCellRangeMoc2Iter<u64, Time<u64>, u64, Hpx::<u64>>) -> Self {
        let (_, upp) = it.size_hint();
        let ub = upp.unwrap_or(100);
        let mut t: Vec::<Range<u64>> = Vec::with_capacity(ub);
        let mut s: Vec::<Ranges<u64>> = Vec::with_capacity(ub);
        for elem in it {
            let (moc_t, moc_s) = elem.mocs();
            /* Simpler but we want to avoid the copy of the s_moc for the last t_range
            for range_t in moc_t.into_cellcellrange_moc_iter().ranges().peekable() {
                t.push(range_t);
                s.push(moc_s.moc_ranges().ranges().clone())
            }*/
            let sranges = Ranges::new_unchecked(moc_s.into_cellcellrange_moc_iter().ranges().collect());
            let mut it = moc_t.into_cellcellrange_moc_iter().ranges().peekable();
            while it.peek().is_some() {
                let range_t = it.next().unwrap();
                t.push(range_t);
                s.push(sranges.clone())
            }
            if let Some(range_t) = it.next() {
                t.push(range_t);
                s.push(sranges)
            }
        }
        HpxRanges2D(Moc2DRanges::<u64, Time<u64>, u64, Hpx<u64>>::new(t, s))
    }
}

impl From<CellMoc2Iter<u64, Time<u64>, u64, Hpx::<u64>>> for HpxRanges2D<u64, Time<u64>, u64> {

    fn from(it: CellMoc2Iter<u64, Time<u64>, u64, Hpx::<u64>>) -> Self {
        let (_, upp) = it.size_hint();
        let ub = upp.unwrap_or(100);
        let mut t: Vec::<Range<u64>> = Vec::with_capacity(ub);
        let mut s: Vec::<Ranges<u64>> = Vec::with_capacity(ub);
        for elem in it {
            let (moc_t, moc_s) = elem.mocs();
            /* Simpler but we want to avoid the copy of the s_moc for the last t_range
            for range_t in moc_t.into_cell_moc_iter().ranges().peekable() {
                t.push(range_t);
                s.push(moc_s.moc_ranges().ranges().clone())
            }*/
            let sranges = Ranges::new_unchecked(moc_s.into_cell_moc_iter().ranges().collect());
            let mut it = moc_t.into_cell_moc_iter().ranges().peekable();
            while it.peek().is_some() {
                let range_t = it.next().unwrap();
                t.push(range_t);
                s.push(sranges.clone())
            }
            if let Some(range_t) = it.next() {
                t.push(range_t);
                s.push(sranges)
            }
        }
        HpxRanges2D(Moc2DRanges::<u64, Time<u64>, u64, Hpx<u64>>::new(t, s))
    }
}




// Adaptor to write FITs
pub struct TimeSpaceRangesIter<'a, T: Idx> {
    depth_max_t: u8,
    depth_max_s: u8,
    it_t: Peekable<slice::Iter<'a, Range<T>>>,
    it_s: Peekable<slice::Iter<'a, Ranges<T>>>
}
impl<'a, T: Idx> HasTwoMaxDepth for TimeSpaceRangesIter<'a, T> {
    fn depth_max_1(&self) -> u8 {
        self.depth_max_t
    }
    fn depth_max_2(&self) -> u8 {
        self.depth_max_s
    }
}
impl<'a, T: Idx> ZSorted for TimeSpaceRangesIter<'a, T> {}
impl<'a, T: Idx> NonOverlapping for TimeSpaceRangesIter<'a, T> {}
impl<'a, T: Idx> MOC2Properties for TimeSpaceRangesIter<'a, T> {}
impl<'a, T: Idx> Iterator for TimeSpaceRangesIter<'a, T> {
    type Item = RangeMOC2Elem<T, Time<T>, T, Hpx<T>>;
    fn next(&mut self) -> Option<Self::Item> {
        if let (Some(t_range), Some(s_ranges)) = (self.it_t.next(), self.it_s.next()) {
            let mut t = Vec::<Range<T>>::new();
            t.push(t_range.clone());
            while let Some(next_s_ranges) = self.it_s.peek() {
                if next_s_ranges == &s_ranges {
                    t.push(self.it_t.next().unwrap().clone());
                    self.it_s.next().unwrap();
                } else {
                    break;
                }
            }
            Some(RangeMOC2Elem::new(
                RangeMOC::new(self.depth_max_t, Ranges::new_unchecked(t).into()),
                RangeMOC::new(self.depth_max_s,s_ranges.clone().into())
            ))
        } else {
            None
        }
    }
}
impl<'a, T: Idx> RangeMOC2Iterator<
    T, Time<T>, RangeMocIter<T, Time<T>>,
    T,  Hpx<T>, RangeMocIter<T, Hpx<T>>,
    RangeMOC2Elem<T, Time<T>, T, Hpx<T>>
> for TimeSpaceRangesIter<'a, T> {}

/*
#[cfg(test)]
mod tests {
    use crate::nestedranges2d::HpxRanges2D;

    use num::{Integer, PrimInt};
    use std::ops::Range;
}*/
