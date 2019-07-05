use std::mem;
use std::ops::Range;

pub fn to_flat_vec<T>(mut data: Vec<Range<T>>) -> Vec<T> {
    let len = data.len() << 1;
    let cap = data.capacity();
    let ptr = data.as_mut_ptr() as *mut T;

    mem::forget(data);

    let result = unsafe {
        Vec::from_raw_parts(ptr, len, cap)
    };

    result
}