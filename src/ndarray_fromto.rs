//! This module contains conversion from/to ndarray.
//! 
//! # Remark
//! We cannot use the [From](https://doc.rust-lang.org/std/convert/trait.From.html) trait here
//! since both ndarray objects and converted objects are not defined in this crates 
//! (and we wanted to remove the ndarray dependency from the MOC crates).

use std::{
  mem,
  ops::Range,
};

use num::One;

use ndarray::{Array1, Array2};

use moc::{
  utils,
  idx::Idx,
  qty::MocQty,
  ranges::{Ranges, SNORanges},
  elemset::range::{
    MocRanges,
    uniq::HpxUniqRanges
  },
  hpxranges2d::HpxRanges2D
};


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

pub fn array2_to_vec_ranges<T: Idx>(mut input: Array2<T>) -> Vec<Range<T>> {
  let shape = input.shape();
  if shape[1] != 2 {
    // Unrecognized shape of a Array2 coming
    // from MOCPy python code
    unreachable!("Unrecognized Array2 shape coming from  MOCPy python code {:?}", shape);
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
