//! This module contains the structures used to perform operations on Range MOC iterators.

pub mod check;
pub mod convert;

// pub mod not;
// pub mod degrade;

pub mod and;   // <=> intersection
pub mod or;    // <=> union
pub mod minus; // <=> mocpy difference = Aladin Soustracction
pub mod xor;   // <=> Aladin Difference