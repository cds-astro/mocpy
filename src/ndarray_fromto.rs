//! This module contains conversion from/to ndarray.
//! 
//! # Remark
//! We cannot use the [From](https://doc.rust-lang.org/std/convert/trait.From.html) trait here
//! since both ndarray objects and converted objects are not defined in this crates 
//! (and we wanted to remove the ndarray dependency from the MOC crates).

use std::mem;
use std::ops::Range;

use num::One;

use ndarray::{Array1, Array2};

use moc::utils;
use moc::idx::Idx;
use moc::qty::{MocQty, Hpx};
use moc::ranges::{Ranges, SNORanges};
use moc::elemset::range::{
  MocRanges,
  uniq::HpxUniqRanges
};
use moc::hpxranges2d::HpxRanges2D;
use moc::mocranges2d::Moc2DRanges;

pub fn uniq_ranges_to_array1<T: Idx>(input: HpxUniqRanges<T>) -> Array1<T> {
  // let ranges = input.ranges;

  let mut result: Vec<T> = Vec::<T>::new();

  // Add the UNIQs contained into the spatial MOC
  for range in input.iter() {
    for uniq_range in num::range_step(range.start, range.end, One::one()) {
      result.push(uniq_range);
    }
  }
  // Get an Array1 from the Vec<i64> without copying any data
  let result: Array1<T> = result.into();
  result.to_owned()
}

/*
impl From<HpxUniqRanges<u64>> for Array1<u64> {
    fn from(input: HpxUniqRanges<u64>) -> Self {
        uniq_ranges_to_array1d(input)
    }
}
impl From<HpxUniqRanges<i64>> for Array1<i64> {
    fn from(input: HpxUniqRanges<i64>) -> Self {
        uniq_ranges_to_array1d(input)
    }
}
*/


pub fn ranges_to_array2<T: Idx>(input: Ranges<T>) -> Array2<T> {
  vec_range_to_array2(input.0.to_vec())
}

pub fn vec_range_to_array2<T: Idx>(mut ranges: Vec<Range<T>>) -> Array2<T> {
  // Cast Vec<Range<u64>> to Vec<u64>
  let len = ranges.len();
  let data = utils::flatten(&mut ranges);

  // Get a Array1 from the Vec<u64> without copying any data
  let result: Array1<T> = data.into();

  // Reshape the result to get a Array2 of shape (N x 2) where N is the number
  // of HEALPix cell contained in the moc
  result.into_shape((len, 2)).unwrap().to_owned()
}

pub fn mocranges_to_array2<T: Idx, Q: MocQty<T>>(input: MocRanges<T, Q>) -> Array2<T> {
  ranges_to_array2(input.0)
}

/*
impl<Q: MocQty<u64>> From<MocRanges<u64, Q>> for Array2<u64> {
  fn from(input: MocRanges<u64, Q>) -> Self {
    ranges_to_array2d(input.0)
  }
}
impl<Q: MocQty<i64>> From<MocRanges<i64, Q>> for Array2<i64> {
  fn from(input: MocRanges<i64, Q>) -> Self {
    ranges_to_array2d(input.0)
  }
}
 */

/*
impl From<Ranges<u64>> for Array2<u64> {
  fn from(input: Ranges<u64>) -> Self {
    ranges_to_array2d(input)
  }
}
impl From<Ranges<i64>> for Array2<i64> {
  fn from(input: Ranges<i64>) -> Self {
    ranges_to_array2d(input)
  }
}
 */

pub fn array2_to_vec_ranges<T: Idx>(mut input: Array2<T>) -> Vec<Range<T>> {
  let shape = input.shape();
  if shape[1] != 2 {
    // Unrecognized shape of a Array2 coming
    // from MOCPy python code
    let msg = format!("Unrecognized Array2 shape coming from  MOCPy python code {:?}", shape);
    unreachable!(msg);
  };
  let len = shape[0];
  let ptr = input.as_mut_ptr() as *mut Range<T>;

  mem::forget(input);

  unsafe { Vec::from_raw_parts(ptr, len, len) }
}

/// Convert a Array2 to a Ranges<T>
///
/// This is useful for converting a set of ranges
/// coming from python into a Rust structure
///
/// # Info
///
/// This method is used whenever an operation must be
/// done to a MOC such as logical operations, degradation
/// max depth computation, etc...
pub fn array2_to_ranges<T: Idx>(input: Array2<T>) ->  Ranges<T> {
  let ranges = array2_to_vec_ranges(input);
  Ranges::<T>::new_unchecked(ranges)
}

/// Convert a Array2 to a MocRanges<T, Q>
///
/// This is useful for converting a set of ranges
/// coming from python into a Rust structure
///
/// # Info
///
/// This method is used whenever an operation must be
/// done to a MOC such as logical operations, degradation
/// max depth computation, etc...
pub fn array2_to_mocranges<T: Idx, Q: MocQty<T>>(input: Array2<T>) ->  MocRanges<T, Q> {
  let ranges = array2_to_vec_ranges(input);
  MocRanges::<T, Q>::new_unchecked(ranges)
}

/// Convert a Array2 to a UniqRanges<T>
///
/// This is useful for converting a set of ranges
/// coming from python into a Rust structure
pub fn array2_to_hpx_uniq_ranges<T: Idx, Q: MocQty<T>>(input: Array2<T>) ->   HpxUniqRanges<T> {
  let ranges = array2_to_vec_ranges(input);
  HpxUniqRanges::<T>::new_unchecked(ranges)
}


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
pub fn hpxranges2d_to_array1_i64<T: MocQty<u64>>(input: &HpxRanges2D<u64, T, u64>) ->  Array1<i64> {
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
pub fn try_array1_i64_to_hpx_ranges2<T:MocQty<u64>>(input: Array1<i64>) -> Result<HpxRanges2D<u64, T, u64>, &'static str> {
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
      return Err("Number of time ranges and spatial coverages do not match.");
    }
    // No need to make it consistent because it comes
    // from python
    Moc2DRanges::<u64, T, u64, Hpx<u64>>::new(t, s)
  };

  Ok(HpxRanges2D(ranges))
}


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
pub fn try_array1_u64_to_hpx_ranges2<T:MocQty<u64>>(input: Array1<u64>) -> Result<HpxRanges2D<u64, T, u64>, &'static str> {
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
      return Err("Number of time ranges and spatial coverages do not match.");
    }
    // No need to make it consistent because it comes
    // from python
    Moc2DRanges::<u64, T, u64, Hpx<u64>>::new(t, s)
  };

  Ok(HpxRanges2D(ranges))
}


#[cfg(test)]
mod tests {
  use ndarray::Array2;
  use moc::ranges::Ranges;
  use crate::ndarray_fromto::{ranges_to_array2, array2_to_vec_ranges};

  #[test]
  fn test_empty_array2_to_vec_ranges() {
    let empty_array = Array2::<u64>::zeros((0, 2));

    let result = array2_to_vec_ranges(empty_array);

    assert_eq!(result, vec![]);
  }
  
  #[test]
  fn empty_ranges_to_array2() {
    let ranges = Ranges::<u64>::new_unchecked(vec![]);

    let result = ranges_to_array2(ranges);
    assert_eq!(result, Array2::<u64>::zeros((0, 2)));
  }
}