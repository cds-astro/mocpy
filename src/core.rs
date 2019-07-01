use num::{Integer, PrimInt, One};

use intervals::ranges::{Ranges, Ranges2D};
use intervals::bounded::Bounded;

use rayon::prelude::*;

fn create_coverage<T, S>(x: Vec<T>, y: Vec<S>, d1: i8, d2: i8) -> Result<Ranges2D<T, S>, &'static str>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {

    let s1 = ((<T>::MAXDEPTH - d1) << 1) as u32;
    let mut off1: T = One::one();
    off1 = off1.unsigned_shl(s1) - One::one();

    let mut m1: T = One::one();
    m1 = m1.checked_mul(&!off1).unwrap();

    let x = x.into_par_iter()
             .map(|r| {
                let a: T = r & m1;
                let b: T = r.checked_add(&off1).unwrap();
                a..b
             })
             .collect::<Vec<_>>();

    let s2 = ((<S>::MAXDEPTH - d2) << 1) as u32;
    let mut off2: S = One::one();
    off2 = off2.unsigned_shl(s2);

    let mut m2: S = One::one();
    m2 = m2.checked_mul(&!off2).unwrap();

    let y = y.into_par_iter()
             .map(|r| {
                 let a: S = r & m2;
                 let b: S = r.checked_add(&off2).unwrap();
                 Ranges::<S>::new(vec![a..b], None, false)
             })
             .collect::<Vec<_>>();

    Ranges2D::<T, S>::new(x, y, None, true)
}

// Returns all the 2nd dim ranges for which their 1st dim range
// is contained into ``x``
//
// This method checks for all the 1st dim ranges that
// lie into the range set ``x``.
// It then returns the union of the 2nd dim set of ranges
// of the matching ranges (i.e. the ones contained into the
// range set ``x``).
fn project_on_second_dim<T, S>(x: &Ranges<T>, coverage: &Ranges2D<T, S>) -> Ranges<S>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
    let result = coverage.x
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
            || Ranges::<S>::new(vec![], None, false),
            |s1, s2| {
                s1.union(&s2)
            }
        );
    result
}

// Returns all the 1st dim ranges for which their 2nd dim set of ranges
// are contained into ``y``
//
// Works in the same manner of ``project_on_first_dim``
fn project_on_first_dim<T, S>(y: &Ranges<S>, coverage: &Ranges2D<T, S>) -> Ranges<T>
where T: Integer + PrimInt + Bounded<T> + Send + Sync + std::fmt::Debug,
      S: Integer + PrimInt + Bounded<S> + Send + Sync + std::fmt::Debug {
    let t_ranges = coverage.x
        .par_iter()
        .zip_eq(coverage.y.par_iter())
        // Filter the time ranges to keep only those
        // that lie into ``x``
        .filter_map(|(t, s)| {
            for r in s.0.iter() {
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

    Ranges::<T>::new(t_ranges, None, true)
}
