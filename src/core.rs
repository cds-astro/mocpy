use ndarray::{Zip, Array, Array1, Array2, Axis};
use ndarray_parallel::prelude::*;
use num::{Integer, PrimInt, One};

use intervals::intervals::Intervals;
use intervals::ranges::{Ranges, Ranges2D};
use intervals::bounded::Bounded;

use std::ops::Range;

fn create_coverage<T, S>(x: Vec<T>, y: Vec<S>, d1: i8, d2: i8) -> Result<Ranges2D<T, S>, &'static str>
where T: Integer + PrimInt + Bounded<T> + Send + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + std::fmt::Debug {

    let shift1 = ((<T>::MAXDEPTH - d1) << 1) as u32;
    let mut offset1: T = One::one();
    offset1.unsigned_shl(shift1);

    let x = x.into_iter()
             .map(|r| {
                 let t1 = r.unsigned_shr(shift1);
                 t1..(t1 + offset1)
             })
             .collect::<Vec<_>>();

    let shift2 = ((<S>::MAXDEPTH - d2) << 1) as u32;
    let mut offset2: S = One::one();
    offset2 = offset2.unsigned_shl(shift2);

    let y = y.into_iter()
             .map(|r| {
                 let pix = r.unsigned_shr(shift2);
                 Ranges::<S>::new(vec![pix..(pix + offset2)], None, false)
             })
             .collect::<Vec<_>>();
    
    Ranges2D::<T, S>::new(x, y, None, true)?
}

fn project_on_first_dim<T, S>(x: &Ranges<T>, coverage: &Ranges2D<T, S>) -> Option<Ranges<S>> {
    coverage.x
	.par_iter()
	.zip(coverage.y.par_iter())
    	// Filter the time ranges to keep only those
	// that lie into ``x``
	.filter(|t, s| {
	    x.contains(t)  
	})
	// Compute the union of all the 2nd dim ranges
	// that have been kept
	.reduce_with(|(t1, ref s1), (t2, ref s2)| {
	    s1.union(s2)
	})
}
