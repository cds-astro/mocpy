
use std::slice;
use std::iter::Peekable;
use std::ops::Range;
use std::convert::From;

use num::One;
use rayon::prelude::*;

use crate::idx::Idx;
use crate::moc::{ZSorted, NonOverlapping, RangeMOCIntoIterator, CellMOCIterator, CellMOCIntoIterator, CellOrCellRangeMOCIterator, CellOrCellRangeMOCIntoIterator, range::{RangeMOC, RangeMocIter}, RangeMOCIterator};
use crate::moc2d::{HasTwoMaxDepth, MOC2Properties, RangeMOC2Iterator, range::RangeMOC2Elem, RangeMOC2ElemIt};
use crate::qty::{MocQty, Hpx, Time};
use crate::ranges::{SNORanges, ranges2d::SNORanges2D, Ranges};
use crate::ranges::ranges2d::Ranges2D;
use crate::elemset::range::{MocRanges, HpxRanges};
use crate::mocranges2d::Moc2DRanges;
use crate::moc2d::cell::CellMoc2Iter;
use crate::moc2d::cellcellrange::CellOrCellRangeMoc2Iter;
/// Declaration of the ST-MOC type
pub type TimeSpaceMoc<T, S> = HpxRanges2D::<T, Time<T>, S>;

// Just to be able to define specific methods on this struct
#[derive(Debug)]
pub struct HpxRanges2D<TT: Idx, T: MocQty<TT>, ST: Idx>(pub Moc2DRanges<TT, T, ST, Hpx<ST>>);

impl<TT: Idx> HpxRanges2D<TT, Time<TT>, TT> {
    pub fn time_space_iter(&self, depth_max_t: u8, depth_max_s: u8) -> TimeSpaceRangesIter<'_, TT> {
        // let (depth_max_t, depth_max_s) = self.compute_min_depth();
        TimeSpaceRangesIter {
            depth_max_t,
            depth_max_s,
            it_t: self.0.ranges2d.x.iter().peekable(),
            it_s: self.0.ranges2d.y.iter().peekable(),
        }
    }
    
   
    pub fn from_ranges_it<I>(it: I) -> Self
        where I: RangeMOC2Iterator<
                TT, Time::<TT>, RangeMocIter<TT, Time::<TT>>,
                TT, Hpx::<TT>, RangeMocIter<TT, Hpx::<TT>>,
                RangeMOC2Elem<TT, Time::<TT>, TT, Hpx::<TT>>
            >
    {
        let mut t = Vec::<Range<TT>>::new();
        let mut s = Vec::<Ranges<TT>>::new();
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
        HpxRanges2D(Moc2DRanges::<TT, Time<TT>, TT, Hpx<TT>>::new(t, s))
    }

    pub fn from_ranges_it_gen<I, J, K, L>(it: L) -> Self
        where
          I: RangeMOCIterator<TT, Qty=Time::<TT>>,
          J: RangeMOCIterator<TT, Qty=Hpx::<TT>>,
          K: RangeMOC2ElemIt<TT, Time::<TT>, TT, Hpx::<TT>, It1=I, It2=J>,
          L: RangeMOC2Iterator<
              TT, Time::<TT>, I,
              TT, Hpx::<TT>, J,
              K
          >
    {
        let mut t = Vec::<Range<TT>>::new();
        let mut s = Vec::<Ranges<TT>>::new();
        for elem in it {
            let (moc_t_it, moc_s_it) = elem.range_mocs_it();
            let moc_t = moc_t_it.into_range_moc();
            let moc_s = moc_s_it.into_range_moc();
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
        HpxRanges2D(Moc2DRanges::<TT, Time<TT>, TT, Hpx<TT>>::new(t, s))
    }
}


impl<TT, T, ST> Default for HpxRanges2D<TT, T, ST>
    where
        TT: Idx,
        T: MocQty<TT>,
        ST: Idx,
{
    /// Create a new empty `NestedRanges2D<T, S>`
    fn default() -> HpxRanges2D<TT, T, ST> {
        let ranges = Moc2DRanges::new(vec![], vec![]);
        HpxRanges2D(ranges)
    }
}

impl<TT, T, ST> HpxRanges2D<TT, T, ST>
where
    TT: Idx,
    T: MocQty<TT>,
    ST: Idx,
{
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
                Ranges::<ST>::default,
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

    /*/// Returns the union of the ranges along the `T` axis for which their
    /// `S` ranges intersect ``y``
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
    pub fn project_on_first_dim_v2(
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
    }*/

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
{ }

// The 3 following From contains code redundancy. We probably should do something!

/*impl<I: RangeMOC2Iterator<
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
}*/
/*
impl<T: Idx, I: RangeMOC2Iterator<
    T, Time::<T>, RangeMocIter<T, Time::<T>>,
    T, Hpx::<T>, RangeMocIter<T, Hpx::<T>>,
    RangeMOC2Elem<T, Time::<T>, T, Hpx::<T>>>
> From<I> for HpxRanges2D<T, Time<T>, T> {

    fn from(it: I) -> Self {
        let mut t = Vec::<Range<T>>::new();
        let mut s = Vec::<Ranges<T>>::new();
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
        HpxRanges2D(Moc2DRanges::<T, Time<T>, T, Hpx<T>>::new(t, s))
    }
}*/

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
            let mut t = vec![t_range.clone()];
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
