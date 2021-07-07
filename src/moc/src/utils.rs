
use std::mem;
use std::ops::Range;

pub fn flatten<T>(input: &mut Vec<Range<T>>) -> Vec<T> {
    let mut owned_input = Vec::<Range<T>>::new();
    // We swap the content referred by input with a new
    // allocated vector.
    // This fix the problem when ``input`` is freed by reaching out
    // the end of the caller scope.
    mem::swap(&mut owned_input, input);

    let len = owned_input.len() << 1;
    let cap = len;
    let ptr = owned_input.as_mut_ptr() as *mut T;

    mem::forget(owned_input);

    unsafe { Vec::from_raw_parts(ptr, len, cap) }
}

pub fn unflatten<T>(input: &mut Vec<T>) -> Vec<Range<T>> {
    let mut owned_input = Vec::<T>::new();
    // We swap the content refered by input with a new
    // allocated vector.
    // This fix the problem when ``input`` is freed by reaching out
    // the end of the caller scope.
    mem::swap(&mut owned_input, input);

    let len = owned_input.len() >> 1;
    let cap = len;
    let ptr = owned_input.as_mut_ptr() as *mut Range<T>;

    mem::forget(owned_input);

    unsafe { Vec::from_raw_parts(ptr, len, cap) }
}

#[cfg(test)]
mod tests {

    use std::ops::Range;
    use crate::utils::{flatten, unflatten};
    
    #[test]
    fn test_empty_flatten() {
        let mut empty_ranges = Vec::<Range<u64>>::new();

        let result = flatten(&mut empty_ranges);

        assert_eq!(result, Vec::<u64>::default());
    }

    #[test]
    fn test_empty_unflatten() {
        let mut empty_vec = Vec::<u64>::new();

        let result = unflatten(&mut empty_vec);

        assert_eq!(result, vec![]);
    }
}
