//! Multi range set operations

// or / and / xor / (difference?) of n sets of ranges a the same time
// using the sweep line algo.
// Allows 8 mocs a the same time, using the 8 u8 bits to tell if the MOC at the given
// bit index intersects or not with the sweep line.
