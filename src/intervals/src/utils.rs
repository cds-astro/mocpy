use std::mem;
use std::ops::Range;

pub fn flatten<T>(input: &mut Vec<Range<T>>) -> Vec<T> {
    let mut owned_input = Vec::<Range<T>>::new();
    // We swap the content refered by input with a new
    // allocated vector.
    // This fix the problem when ``input`` is freed by reaching out
    // the end of the caller scope.
    std::mem::swap(&mut owned_input, input);

    let len = owned_input.len() << 1;
    let cap = owned_input.capacity();
    let ptr = owned_input.as_mut_ptr() as *mut T;

    mem::forget(owned_input);

    let result = unsafe { Vec::from_raw_parts(ptr, len, cap) };

    result
}

pub fn unflatten<T>(input: &mut Vec<T>) -> Vec<Range<T>> {
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

    let result = unsafe { Vec::from_raw_parts(ptr, len, cap) };

    result
}

use crate::bounded::Bounded;
use ndarray::Array2;
use num::{Integer, PrimInt};

pub fn array2_to_vec_ranges<T>(mut input: Array2<T>) -> Vec<Range<T>>
where
    T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
{
    let shape = input.shape();
    // Warning: empty Array2 objects coming from MOCPy
    // do have a shape equal to (1, 0).
    let len = if shape == [1, 0] {
        // It is empty
        0
    } else if shape[1] == 2 {
        shape[0]
    } else {
        // Unrecognized shape of a Array2 coming
        // from MOCPy python code
        let msg = format!(
            "Unrecognized Array2 shape coming from \
             MOCPy python code {:?}",
            shape
        );
        unreachable!(msg);
    };
    let cap = len;
    let ptr = input.as_mut_ptr() as *mut Range<T>;

    mem::forget(input);

    unsafe { Vec::from_raw_parts(ptr, len, cap) }
}

#[cfg(test)]
mod tests {
    use crate::utils::array2_to_vec_ranges;
    use ndarray::Array2;
    #[test]
    fn test_empty_array2_to_vec_ranges() {
        let empty_array = Array2::<u64>::zeros((1, 0));

        let result = array2_to_vec_ranges(empty_array);

        assert_eq!(result, vec![]);
    }

    use crate::utils::flatten;
    use std::ops::Range;
    #[test]
    fn test_empty_flatten() {
        let mut empty_ranges = Vec::<Range<u64>>::new();

        let result = flatten(&mut empty_ranges);

        assert_eq!(result, vec![]);
    }

    use crate::utils::unflatten;
    #[test]
    fn test_empty_unflatten() {
        let mut empty_vec = Vec::<u64>::new();

        let result = unflatten(&mut empty_vec);

        assert_eq!(result, vec![]);
    }
}
