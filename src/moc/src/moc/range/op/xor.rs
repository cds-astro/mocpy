
use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};
use std::cmp::Ordering;

/// Performs a logical `OR` between the two input iterators of ranges.
pub fn xor<T, Q, I1, I2>(
    left_it: I1,
    right_it: I2
) -> XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{
    XorRangeIter::new(left_it, right_it)
}

/// Performs an `OR` operation between two iterators of ranges on-the-fly, while iterating.
pub struct XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{
    last: Option<Range<T>>,
    left_it: I1,
    right_it: I2,
    left: Option<Range<T>>,
    right: Option<Range<T>>,
}

impl <T, Q, I1, I2> XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{

    fn new(mut left_it: I1, mut right_it: I2) -> XorRangeIter<T, Q, I1, I2> {
        let last = if let (Some(range1), Some(range2)) = (left_it.peek_last(), right_it.peek_last()) {
            if range2.end < range1.start { // |--2--| |--1--|
                Some(range1.clone())
            } else if range1.end < range2.start { // |--1--| |--2--|
                Some(range2.clone())
            } else { // Overlapping ranges
                let start = range1.start.min(range2.start);
                let end = range1.end.max(range2.end);
                Some(start..end)
            }
        } else {
            None
        };
        let left = left_it.next();
        let right = right_it.next();
        // Normal behaviour
        XorRangeIter {
            last,
            left_it,
            right_it,
            left,
            right,
        }
    }
}

impl<T, Q, I1, I2> HasMaxDepth for  XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{
    fn depth_max(&self) -> u8 {
        u8::max(self.left_it.depth_max(), self.right_it.depth_max())
    }
}
impl<T, Q, I1, I2> ZSorted for  XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{ }
impl<T, Q, I1, I2> NonOverlapping for XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{ }
impl<T, Q, I1, I2> MOCProperties for XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{ }
impl<T, Q, I1, I2> Iterator for XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{

    type Item = Range<T>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // We can probably simplify this, check again another day!!
            match (&mut self.left, &mut self.right) {
                (Some(ref mut l), Some(ref mut r)) =>
                    if l.end == r.start {        // L--LR--R => return R----R; next L
                        r.start = l.start;
                        self.left = self.left_it.next();
                    } else if r.end == l.start { // R--RL--L => return L----L; next R
                        l.start = r.start;
                        self.right = self.right_it.next();
                    } else if l.end < r.start {  // L--L R--R => return L--L; next L
                        return std::mem::replace(&mut self.left, self.left_it.next());
                    } else if r.end < l.start {  // R--R L--L => return R--R, and next R
                        return std::mem::replace(&mut self.right, self.right_it.next());
                    } else if l.end == r.end {
                        match l.start.cmp(&r.start) {
                            Ordering::Equal => { // LR--LR => next L and next R
                                self.left = self.left_it.next();
                                self.right = self.right_it.next();
                            },
                            Ordering::Less => { // L--R--LR => return L--R; next L next R
                                l.end = r.start;
                                self.right = self.right_it.next();
                                return std::mem::replace(&mut self.left, self.left_it.next());
                            },
                            Ordering::Greater => { // R--L--LR =>return R--L; next L next R
                                r.end = l.start;
                                self.left = self.left_it.next();
                                return std::mem::replace(&mut self.right, self.right_it.next());
                            },
                        }
                    } else if l.end < r.end { //    RL--L--R or R--L--L--R or L--R--L--R
                        match l.start.cmp(&r.start) {
                            Ordering::Equal => { // RL--L--R, R--R, next L
                                r.start = l.end;
                                self.left = self.left_it.next();
                            },
                            Ordering::Less => { // L--R--L--R, L--LR--R return L increment L
                                std::mem::swap(&mut l.end, &mut r.start);
                                return std::mem::replace(&mut self.left, self.left_it.next());
                            },
                            Ordering::Greater => { // R--L--L--R, L--L R--R  return L increment L
                                std::mem::swap(&mut l.end, &mut r.start);
                                std::mem::swap(&mut l.end, &mut l.start);
                                return std::mem::replace(&mut self.left, self.left_it.next());
                            },
                        }
                    } else {
                        assert!(r.end < l.end); //    LR--R--L or L--R--R--L or R--L--R--L
                        match l.start.cmp(&r.start) {
                            Ordering::Equal => { // LR--R--L, L--L, next R
                                l.start = r.end;
                                self.right = self.right_it.next();
                            },
                            Ordering::Less => { //  L--R--R--L, R--R L--L return R increment R
                                std::mem::swap(&mut r.end, &mut l.start);
                                std::mem::swap(&mut r.end, &mut r.start);
                                return std::mem::replace(&mut self.right, self.right_it.next());
                            },
                            Ordering::Greater => { // R--L--R--L, R--RL--L  return R increment R
                                std::mem::swap(&mut r.end, &mut l.start);
                                return std::mem::replace(&mut self.right, self.right_it.next());
                            },
                        }
                    },
                (Some(_), None) => return std::mem::replace(&mut self.left, self.left_it.next()),
                (None, Some(_)) => return std::mem::replace(&mut self.right, self.right_it.next()),
                (None, None) => return None,
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size_hint_l = self.left_it.size_hint();
        let size_hint_r = self.right_it.size_hint();
        if let (Some(n1), Some(n2)) = (size_hint_l.1, size_hint_r.1) {
            // Worse case:
            // * 1 range of T1 contains (n2-1) T2 ranges
            // * and 1 range of T2 contains the remaining (n1-1) cells of T1
            (0, Some(1 + n1 + n2))
        } else {
            (0, None)
        }
    }

}


impl<T, Q, I1, I2> RangeMOCIterator<T> for XorRangeIter<T, Q, I1, I2>
    where
        T: Idx,
        Q: MocQty<T>,
        I1: RangeMOCIterator<T, Qty=Q>,
        I2: RangeMOCIterator<T, Qty=Q>
{
    type Qty = Q;

    fn peek_last(&self) -> Option<&Range<T>> {
        self.last.as_ref()
    }
}


#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::ops::Range;
    use std::path::PathBuf;
    use std::io::BufReader;

    use crate::qty::Hpx;
    use crate::elem::cell::Cell;
    use crate::moc::range::RangeMOC;
    use crate::deser::fits::{
        from_fits_ivoa,
        MocIdxType, MocQtyType, MocType
    };
    use crate::moc::{
        HasMaxDepth,
        RangeMOCIterator, RangeMOCIntoIterator,
        CellMOCIterator, CellMOCIntoIterator,
        range::op::xor::xor
    };

    fn load_moc(filename: &str) -> RangeMOC<u32, Hpx<u32>> {
        let path_buf1 = PathBuf::from(format!("resources/{}", filename));
        let path_buf2 = PathBuf::from(format!("../resources/{}", filename));
        let file = File::open(&path_buf1).or_else(|_| File::open(&path_buf2)).unwrap();
        let reader = BufReader::new(file);
        match from_fits_ivoa(reader) {
            Ok(MocIdxType::U32(MocQtyType::Hpx(MocType::Ranges(moc)))) => {
                let moc = RangeMOC::new(moc.depth_max(), moc.collect());
                moc
            },
            Ok(MocIdxType::U32(MocQtyType::Hpx(MocType::Cells(moc)))) => {
                let moc = RangeMOC::new(
                    moc.depth_max(),
                    moc.into_cell_moc_iter().ranges().collect()
                );
                moc
            },
            _ => unreachable!(false),
        }
    }

    fn load_mocs() -> (RangeMOC<u32, Hpx<u32>>, RangeMOC<u32, Hpx<u32>>) {
        let sdss = load_moc("V_147_sdss12.moc.fits");
        let other = load_moc("CDS-I-125A-catalog_MOC.fits");
        (sdss, other)
    }

    /* Not implemented in MOCPy
    fn xor_ranges(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
        let depth = u8::max(moc_l.depth_max(), moc_r.depth_max());
        let ranges_l = moc_l.into_moc_ranges();
        let ranges_r = moc_r.into_moc_ranges();
        RangeMOC::new(depth, ranges_l.difference(&ranges_r))
    }*/

    // we could also perform the operation without having first collected the iterator we obtain from
    // the FITS file
    fn xor_ranges_it(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
        let xor = xor(moc_l.into_range_moc_iter(), moc_r.into_range_moc_iter());
        RangeMOC::new(xor.depth_max(), xor.collect())
    }

    #[test]
    fn test_xor() {
        let (sdss, other) = load_mocs();
        // let xor_1 = xor_ranges(sdss.clone(), other.clone());
        let xor_2 = xor_ranges_it(sdss.clone(), other.clone());
        // println!("XOR 1 range size: {}", xor_1.moc_ranges().0.0.len());
        println!("XOR 2 range size: {}", xor_2.moc_ranges().0.0.len());
        // println!("XOR 1 cell size: {}", (&xor_1).into_range_moc_iter().cells().collect::<Vec<Cell<u32>>>().len());
        println!("XOR 2 cell size: {}", (&xor_2).into_range_moc_iter().cells().collect::<Vec<Cell<u32>>>().len());
        // Comparison made with Aladin
        assert_eq!((&xor_2).into_range_moc_iter().cells().collect::<Vec<Cell<u32>>>().len(), 462928);
        // Ensures the number of ranges is ok
        assert_eq!((&xor_2).into_range_moc_iter().cells().ranges().collect::<Vec<Range<u32>>>().len(), xor_2.moc_ranges().0.0.len());
        // assert_eq!(xor_1.moc_ranges().0.0.len(), xor_2.moc_ranges().0.0.len());
    }
}
