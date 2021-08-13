//! This module contains the structures used to perform operations on Range MOC iterators.

pub mod check;
pub mod convert;
pub mod merge;

pub mod not;   // <=> complement
pub mod degrade;

pub mod and;   // <=> intersection
pub mod or;    // <=> union
pub mod minus; // <=> mocpy difference = Aladin Soustracction
pub mod xor;   // <=> Aladin Difference