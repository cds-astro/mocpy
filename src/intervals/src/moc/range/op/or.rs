
use std::ops::Range;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};

/// Performs a logical `OR` between the two input iterators of ranges.
pub fn or<T, Q, I1, I2>(
  left_it: I1,
  right_it: I2
) -> OrRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{
  OrRangeIter::new(left_it, right_it)
}

/// Performs an `OR` operation between two iterators of ranges on-the-fly, while iterating.
pub struct OrRangeIter<T, Q, I1, I2>
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

impl <T, Q, I1, I2> OrRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{

  fn new(mut left_it: I1, mut right_it: I2) -> OrRangeIter<T, Q, I1, I2> {
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
    OrRangeIter {
      last,
      left_it,
      right_it,
      left,
      right,
    }
  }
}

impl<T, Q, I1, I2> HasMaxDepth for  OrRangeIter<T, Q, I1, I2>
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
impl<T, Q, I1, I2> ZSorted for  OrRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{ }
impl<T, Q, I1, I2> NonOverlapping for OrRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{ }
impl<T, Q, I1, I2> MOCProperties for OrRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{ }
impl<T, Q, I1, I2> Iterator for OrRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{

  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    loop {
      match (&mut self.left, &mut self.right) {
        (Some(ref mut l), Some(ref mut r)) =>
          if l.end < r.start {        // L--L R--R
            return std::mem::replace(&mut self.left, self.left_it.next());
          } else if r.end < l.start { // R--R L--L
            return std::mem::replace(&mut self.right, self.right_it.next());
          } else if l.end <= r.end {  //    R--L--L--R
            if l.start < r.start {    // or L--R--L--R
              r.start = l.start;
            }
            self.left = consume_while_end_lower_than(&mut self.left_it, r.end);
          } else {                    //    L--R--R--L
            if r.start < l.start {    // or R--L--R--L
              l.start = r.start;
            }
            self.right = consume_while_end_lower_than(&mut self.right_it, l.end);
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


impl<T, Q, I1, I2> RangeMOCIterator<T> for OrRangeIter<T, Q, I1, I2>
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

fn consume_while_end_lower_than<T, Q, I>(it: &mut I, to: T) -> Option<Range<T>>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
{
  let mut curr = it.next();
  while let Some(c) = &curr {
    if c.end > to {
      break;
    } else {
      curr = it.next();
    }
  }
  curr
}




#[cfg(test)]
mod tests {
  use std::fs::File;
  use std::io::BufReader;
  use std::path::PathBuf;

  use crate::qty::Hpx;
  use crate::moc::range::RangeMOC;
  use crate::deser::fits::{
    from_fits_ivoa,
    MocIdxType, MocQtyType, MocType
  };
  use crate::moc::{HasMaxDepth, RangeMOCIntoIterator, CellMOCIntoIterator, CellMOCIterator};
  use crate::moc::range::op::or::or;
  use crate::ranges::SNORanges;

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

  fn or_ranges(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
    let depth = u8::max(moc_l.depth_max(), moc_r.depth_max());
    let ranges_l = moc_l.into_moc_ranges();
    let ranges_r = moc_r.into_moc_ranges();
    RangeMOC::new(depth, ranges_l.union(&ranges_r))
  }

  // we could also perform the operation without having first collected the iteartor we obtain from
  // the FITS file
  fn or_ranges_it(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
    let or = or(moc_l.into_range_moc_iter(), moc_r.into_range_moc_iter());
    RangeMOC::new(or.depth_max(), or.collect())
  }

  #[test]
  fn test_or() {
    let (sdss, other) = load_mocs();
    let and_1 = or_ranges(sdss.clone(), other.clone());
    let and_2 = or_ranges_it(sdss.clone(), other.clone());
    println!("OR size: {}", and_1.moc_ranges().0.0.len());
    assert_eq!(and_1.moc_ranges().0.0.len(), and_2.moc_ranges().0.0.len());
  }
}
