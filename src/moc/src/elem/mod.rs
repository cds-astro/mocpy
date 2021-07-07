
// cell
// ucell
// zcell

// range
// cellrange

// zcell or range
// - 10 -> range lower bound
// - 11 -> range upper bound
// - 0  -> zcell (like a range, but << 1 and with a sentinel bit coding the depth)

pub mod cell;
pub mod cellrange;
pub mod cellcellrange;
pub mod range;
pub mod valuedcell;
