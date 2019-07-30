use std::mem;
use std::ops::Range;

pub fn flatten<T>(mut data: Vec<Range<T>>) -> Vec<T> {
    let len = data.len() << 1;
    let cap = data.capacity();
    let ptr = data.as_mut_ptr() as *mut T;

    mem::forget(data);

    let result = unsafe {
        Vec::from_raw_parts(ptr, len, cap)
    };

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

    let result = unsafe {
        Vec::from_raw_parts(ptr, len, cap)
    };
    
    result
}